#define _GNU_SOURCE
#include <errno.h>
#include <linux/perf_event.h>
#include <stdio.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
int main(void){const char*names[]={"cpu_cycles","instructions","cache_misses","branch_misses"};const unsigned long long configs[]={PERF_COUNT_HW_CPU_CYCLES,PERF_COUNT_HW_INSTRUCTIONS,PERF_COUNT_HW_CACHE_MISSES,PERF_COUNT_HW_BRANCH_MISSES};puts("{\"scope\":\"own-process user-space PMU capability probe; not campaign timing\",\"events\":[");for(unsigned i=0;i<4;i++){struct perf_event_attr a;memset(&a,0,sizeof a);a.type=PERF_TYPE_HARDWARE;a.size=sizeof a;a.config=configs[i];a.disabled=1;a.exclude_kernel=1;a.exclude_hv=1;int fd=(int)syscall(__NR_perf_event_open,&a,0,-1,-1,0);if(fd<0){printf("{\"name\":\"%s\",\"available\":false,\"value\":null,\"errno\":%d,\"reason\":\"%s\"}%s\n",names[i],errno,strerror(errno),i==3?"":",");continue;}unsigned long long value=0;int ok=ioctl(fd,PERF_EVENT_IOC_RESET,0)==0&&ioctl(fd,PERF_EVENT_IOC_ENABLE,0)==0;volatile unsigned long long sum=0;for(unsigned k=0;k<1000000;k++)sum+=k;ok=ok&&ioctl(fd,PERF_EVENT_IOC_DISABLE,0)==0&&read(fd,&value,sizeof value)==sizeof value;close(fd);printf("{\"name\":\"%s\",\"available\":%s,\"probe_value\":%llu,\"not_a_solver_measurement\":true}%s\n",names[i],ok?"true":"false",value,i==3?"":",");}puts("]}");return 0;}
