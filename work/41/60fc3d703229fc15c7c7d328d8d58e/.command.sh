#!/bin/bash -ue
samtools         view         -@64         test.bam         | grep -m 1 MM:Z
