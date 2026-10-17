#!/usr/bin/env python3
"""Run the Lambda broker and validator as one small external control service."""

import hashlib
import json
import os
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import boto3
from botocore.exceptions import ClientError

from . import broker, validator


PORT = int(os.environ.get("PORT", "8080"))
QUEUE_URL = os.environ["VALIDATOR_QUEUE_URL"]
MAX_REQUEST_BYTES = 1024 * 1024
RATE_PER_MINUTE = int(os.environ.get("BROKER_RATE_PER_MINUTE", "300"))
hits = {}
hitsLock = threading.Lock()


def rateAllowed(header):
    digest = hashlib.sha256((header or "").encode()).hexdigest()
    minute = int(time.time()) // 60
    with hitsLock:
        oldMinute, count = hits.get(digest, (minute, 0))
        if oldMinute != minute:
            oldMinute, count = minute, 0
        count += 1
        hits[digest] = (oldMinute, count)
        if len(hits) > 10000:
            stale = [key for key, value in hits.items() if value[0] < minute - 2]
            for key in stale:
                hits.pop(key, None)
        return count <= RATE_PER_MINUTE


class Handler(BaseHTTPRequestHandler):
    server_version = "ecc2k130-upload-gateway/1"

    def log_message(self, template, *args):
        print("http " + template % args, flush=True)

    def sendJson(self, status, body):
        raw = json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
        self.send_response(status)
        self.send_header("content-type", "application/json")
        self.send_header("cache-control", "no-store")
        self.send_header("content-length", str(len(raw)))
        self.end_headers()
        self.wfile.write(raw)

    def do_GET(self):
        if self.path == "/health":
            self.sendJson(200, {"ok": True})
        else:
            self.sendJson(404, {"error": "route not found"})

    def do_POST(self):
        if not rateAllowed(self.headers.get("authorization")):
            self.sendJson(429, {"error": "worker rate limit exceeded"})
            return
        try:
            length = int(self.headers.get("content-length", "0"))
        except ValueError:
            length = -1
        if length < 0 or length > MAX_REQUEST_BYTES:
            self.sendJson(413, {"error": "request body is too large"})
            return
        raw = self.rfile.read(length)
        event = {"rawPath": self.path.split("?", 1)[0],
                 "requestContext": {"http": {"method": "POST"}},
                 "headers": {key.lower(): value for key, value in self.headers.items()},
                 "body": raw.decode(errors="strict"), "isBase64Encoded": False}
        result = broker.lambdaHandler(event, None)
        self.sendJson(int(result["statusCode"]), json.loads(result["body"]))


def validationLoop():
    queue = boto3.client("sqs", region_name=os.environ.get("AWS_REGION", "us-west-2"))
    while True:
        try:
            page = queue.receive_message(QueueUrl=QUEUE_URL, MaxNumberOfMessages=1,
                                         WaitTimeSeconds=20, VisibilityTimeout=960)
            for message in page.get("Messages", []):
                validator.lambdaHandler({"Records": [{"body": message["Body"]}]}, None)
                queue.delete_message(QueueUrl=QUEUE_URL, ReceiptHandle=message["ReceiptHandle"])
        except Exception as error:
            detail = ""
            if isinstance(error, ClientError):
                detail = ": " + error.response.get("Error", {}).get("Code", "AWS error")
            print("validator poll failed: %s%s" % (type(error).__name__, detail), flush=True)
            time.sleep(5)


def main():
    thread = threading.Thread(target=validationLoop, daemon=True)
    thread.start()
    server = ThreadingHTTPServer(("0.0.0.0", PORT), Handler)
    print("broker listening on port %d" % PORT, flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
