# fastx-mutate-tools

In this jar I have included a few simple tools that permit to introduce SNPs/indels/bisulfite-conversions in fasta/fastq files. Be aware that mutations are inserted uniformly according to the specified probabilities, so the tools do not simulate realistic mutations. However, you may find them useful to test performances of DNA/bisulfite aligners. The tools are designed to be coupled with read simulators such as, e.g., wgsim or SimSeq.

Usage: execute as

> java -jar FastxMutateTools.jar

for info usage
