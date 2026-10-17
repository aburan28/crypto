#include "json.hpp"
#include <chrono>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
using json=nlohmann::json;
int main(int argc,char**argv){
    if(argc<5){std::cerr<<"RECEIPT.json STDOUT STDERR EXE [ARG...]\n";return 2;}
    int receipt=open(argv[1],O_WRONLY|O_CREAT|O_EXCL,0600);
    if(receipt<0){std::cerr<<"refuse receipt overwrite: "<<strerror(errno)<<'\n';return 2;}
    auto start=std::chrono::steady_clock::now();pid_t child=fork();
    if(child==0){
        int out=open(argv[2],O_WRONLY|O_CREAT|O_EXCL,0600),err=open(argv[3],O_WRONLY|O_CREAT|O_EXCL,0600);
        if(out<0||err<0)_exit(125);
        if(dup2(out,STDOUT_FILENO)<0||dup2(err,STDERR_FILENO)<0)_exit(125);
        close(out);close(err);close(receipt);execvp(argv[4],argv+4);perror("execvp");_exit(127);
    }
    if(child<0){close(receipt);return 2;}
    int status=0;rusage ru{};while(wait4(child,&status,0,&ru)<0){if(errno!=EINTR){close(receipt);return 2;}}
    auto elapsed=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-start).count();
    json args=json::array();for(int i=4;i<argc;i++)args.push_back(argv[i]);
    int code=WIFEXITED(status)?WEXITSTATUS(status):128+WTERMSIG(status);
    json record={{"argv",args},{"exit_status",code},{"wall_ns",elapsed},{"child_user_us",ru.ru_utime.tv_sec*1000000LL+ru.ru_utime.tv_usec},{"child_system_us",ru.ru_stime.tv_sec*1000000LL+ru.ru_stime.tv_usec},{"max_rss_kib",ru.ru_maxrss},{"minor_faults",ru.ru_minflt},{"major_faults",ru.ru_majflt},{"voluntary_switches",ru.ru_nvcsw},{"involuntary_switches",ru.ru_nivcsw},{"stdout",argv[2]},{"stderr",argv[3]}};
    std::string data=record.dump(2)+"\n";size_t done=0;while(done<data.size()){ssize_t n=write(receipt,data.data()+done,data.size()-done);if(n<0&&errno==EINTR)continue;if(n<=0){close(receipt);return 2;}done+=size_t(n);}close(receipt);std::cout<<record.dump()<<'\n';return code;
}
