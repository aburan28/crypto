//! AWS clients and the few calls more than one command needs.

use anyhow::{Context, Result};
use aws_config::{BehaviorVersion, Region};
use aws_sdk_ec2::types::Filter;
use aws_smithy_http_client::{Builder, Connector, proxy::ProxyConfig, tls};
use serde_json::Value;
use std::path::Path;

pub struct Cloud {
    pub region: String,
    pub ec2: aws_sdk_ec2::Client,
    pub s3: aws_sdk_s3::Client,
    pub ssm: aws_sdk_ssm::Client,
}

/// EPHEMERAL_REGION, then AWS_REGION, then AWS_DEFAULT_REGION, then us-east-1.
pub fn region() -> String {
    ["EPHEMERAL_REGION", "AWS_REGION", "AWS_DEFAULT_REGION"]
        .iter()
        .find_map(|k| std::env::var(k).ok().filter(|v| !v.is_empty()))
        .unwrap_or_else(|| "us-east-1".into())
}

/// An HTTPS client that honours HTTPS_PROXY/NO_PROXY and AWS_CA_BUNDLE, as the
/// aws CLI does; the SDK's default client does neither.
fn http_client() -> Result<aws_smithy_runtime_api::client::http::SharedHttpClient> {
    let mut trust = tls::TrustStore::default();
    if let Ok(path) = std::env::var("AWS_CA_BUNDLE") {
        let pem = std::fs::read(&path).with_context(|| format!("reading AWS_CA_BUNDLE {path}"))?;
        trust = trust.with_pem_certificate(pem);
    }
    let ctx = tls::TlsContext::builder().with_trust_store(trust).build()?;
    Ok(
        Builder::new().build_with_connector_fn(move |settings, components| {
            let mut b = Connector::builder().proxy_config(ProxyConfig::from_env());
            if let Some(s) = settings {
                b = b.connector_settings(s.clone());
            }
            if let Some(sleep) = components.and_then(|c| c.sleep_impl()) {
                b = b.sleep_impl(sleep);
            }
            b.tls_provider(tls::Provider::Rustls(
                tls::rustls_provider::CryptoMode::AwsLc,
            ))
            .tls_context(ctx.clone())
            .build()
        }),
    )
}

impl Cloud {
    pub async fn new() -> Result<Cloud> {
        let region = region();
        let conf = aws_config::defaults(BehaviorVersion::latest())
            .region(Region::new(region.clone()))
            .http_client(http_client()?)
            .load()
            .await;
        Ok(Cloud {
            region,
            ec2: aws_sdk_ec2::Client::new(&conf),
            s3: aws_sdk_s3::Client::new(&conf),
            ssm: aws_sdk_ssm::Client::new(&conf),
        })
    }

    pub async fn parameter(&self, name: &str) -> Result<String> {
        let out = self
            .ssm
            .get_parameter()
            .name(name)
            .send()
            .await
            .map_err(err)?;
        out.parameter()
            .and_then(|p| p.value())
            .map(String::from)
            .context("parameter has no value")
    }

    pub async fn put(&self, bucket: &str, key: &str, body: Vec<u8>) -> Result<()> {
        self.s3
            .put_object()
            .bucket(bucket)
            .key(key)
            .body(body.into())
            .send()
            .await
            .map_err(err)
            .with_context(|| format!("uploading s3://{bucket}/{key}"))?;
        Ok(())
    }

    /// The object's bytes, or None when it is missing or unreadable.
    pub async fn get(&self, bucket: &str, key: &str, range: Option<String>) -> Option<Vec<u8>> {
        let out = self
            .s3
            .get_object()
            .bucket(bucket)
            .key(key)
            .set_range(range)
            .send()
            .await
            .ok()?;
        Some(out.body.collect().await.ok()?.into_bytes().to_vec())
    }

    pub async fn status(&self, bucket: &str, id: &str) -> Option<Value> {
        let raw = self
            .get(bucket, &format!("jobs/{id}/status.json"), None)
            .await?;
        serde_json::from_slice(&raw).ok()
    }

    /// The newest instance tagged for this job: (instance id, state).
    pub async fn instance_of(&self, id: &str) -> Option<(String, String)> {
        let out = self
            .ec2
            .describe_instances()
            .filters(
                Filter::builder()
                    .name("tag:ephemeral:job")
                    .values(id)
                    .build(),
            )
            .send()
            .await
            .ok()?;
        let mut all: Vec<_> = out
            .reservations()
            .iter()
            .flat_map(|r| r.instances())
            .collect();
        all.sort_by_key(|i| i.launch_time().map(|t| t.secs()).unwrap_or(0));
        let i = all.last()?;
        Some((
            i.instance_id()?.to_string(),
            i.state()
                .and_then(|s| s.name())
                .map(|n| n.as_str().to_string())
                .unwrap_or_default(),
        ))
    }

    /// `aws s3 sync s3://bucket/prefix/ dir/ --exclude runner.sh`: download what is
    /// missing, a different size, or newer than the local copy.
    pub async fn sync_down(&self, bucket: &str, prefix: &str, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;
        let mut pages = self
            .s3
            .list_objects_v2()
            .bucket(bucket)
            .prefix(prefix)
            .into_paginator()
            .send();
        while let Some(page) = pages.next().await {
            for obj in page.map_err(err)?.contents() {
                let Some(key) = obj.key() else { continue };
                let rel = &key[prefix.len()..];
                if rel.is_empty() || rel == "runner.sh" || rel.split('/').any(|p| p == "..") {
                    continue;
                }
                let dest = dir.join(rel);
                if let Ok(meta) = std::fs::metadata(&dest) {
                    let same_size = meta.len() as i64 == obj.size().unwrap_or(-1);
                    let local = meta
                        .modified()
                        .ok()
                        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok());
                    let remote = obj.last_modified().map(|t| t.secs());
                    if same_size
                        && matches!((local, remote), (Some(l), Some(r)) if l.as_secs() as i64 >= r)
                    {
                        continue;
                    }
                }
                let Some(bytes) = self.get(bucket, key, None).await else {
                    eprintln!("ephemeral: could not download s3://{bucket}/{key}");
                    continue;
                };
                if let Some(parent) = dest.parent() {
                    std::fs::create_dir_all(parent)?;
                }
                std::fs::write(&dest, bytes)?;
            }
        }
        Ok(())
    }
}

/// An SDK error with its service code and message, not just "service error".
pub fn err<E: std::error::Error + Send + Sync + 'static>(e: E) -> anyhow::Error {
    anyhow::anyhow!("{}", aws_sdk_ec2::error::DisplayErrorContext(&e))
}
