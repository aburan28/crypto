{{- define "taskq.name" -}}{{ .Release.Name }}-taskq{{- end -}}
{{- define "taskq.labels" -}}
app.kubernetes.io/name: taskq
app.kubernetes.io/instance: {{ .Release.Name }}
app.kubernetes.io/managed-by: {{ .Release.Service }}
helm.sh/chart: {{ .Chart.Name }}-{{ .Chart.Version }}
{{- end -}}
{{- define "taskq.redisSecret" -}}
{{- .Values.redis.existingSecret | default (printf "%s-redis" (include "taskq.name" .)) -}}
{{- end -}}
{{- define "taskq.resultsClaim" -}}
{{- .Values.results.existingClaim | default (printf "%s-results" (include "taskq.name" .)) -}}
{{- end -}}
{{/* Env giving TASKQ_REDIS_URL, from the in-chart Redis or an external URL secret. */}}
{{- define "taskq.redisEnv" -}}
{{- if .Values.redis.enabled }}
- name: REDIS_PASSWORD
  valueFrom: {secretKeyRef: {name: {{ include "taskq.redisSecret" . }}, key: redis-password}}
- name: TASKQ_REDIS_URL
  value: "redis://:$(REDIS_PASSWORD)@{{ include "taskq.name" . }}-redis:6379/0"
{{- else }}
- name: TASKQ_REDIS_URL
  valueFrom: {secretKeyRef: {name: {{ required "redis.externalUrlSecret.name is required when redis.enabled=false" .Values.redis.externalUrlSecret.name }}, key: {{ .Values.redis.externalUrlSecret.key }}}}
{{- end }}
{{- with .Values.namespace }}
- {name: TASKQ_NAMESPACE, value: {{ . | quote }}}
{{- end }}
{{- end -}}
