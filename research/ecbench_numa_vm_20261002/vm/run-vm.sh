#!/bin/sh
# Boot the two-node guest this evidence came from.  Thin orchestration only:
# everything measured runs inside the guest as /init and ecbench.  INIT picks
# the init script: init-v2 (default; with the interruption test) or init.
#   $1: a static x86-64 ecbench binary (see ../README.md for the build)
#   $2: the Ubuntu generic kernel, noble-server-cloudimg-amd64-vmlinuz-generic
#   $3: busybox 1.35.0 x86_64-linux-musl
set -eu
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d "${TMPDIR:-/tmp}/ecbench-vm.XXXXXX")
mkdir -p "$work/root/bin"
cp "$1" "$work/root/ecbench"; cp "$3" "$work/root/bin/busybox"
cp "$here/${INIT:-init-v2}" "$work/root/init"; cp "$here/spec.json" "$work/root/spec.json"
cp "$here/spec-long.json" "$work/root/spec-long.json"
chmod 755 "$work/root/ecbench" "$work/root/bin/busybox" "$work/root/init"
(cd "$work/root" && find . | cpio -o -H newc 2>/dev/null | gzip -9 > "$work/initramfs.cpio.gz")
qemu-system-x86_64 -M q35 -accel tcg,thread=multi -cpu max,vendor=GenuineIntel \
  -smp 8,sockets=2,cores=2,threads=2 -m 1024 \
  -object memory-backend-ram,id=m0,size=512M -object memory-backend-ram,id=m1,size=512M \
  -numa node,nodeid=0,cpus=0-3,memdev=m0 -numa node,nodeid=1,cpus=4-7,memdev=m1 \
  -kernel "$2" -initrd "$work/initramfs.cpio.gz" \
  -append "console=ttyS0 rdinit=/init psi=1 isolcpus=6,7 panic=-1 quiet" \
  -display none -serial "file:$work/console.log" -monitor none -no-reboot -net none
echo "console: $work/console.log"
