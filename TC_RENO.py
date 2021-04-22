import sys, argparse
import subprocess


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Isoform identification')
    required = parser.add_argument_group('required named arguments')
    parser.add_argument('-i',required=True, help='bed files of reads')
    required.add_argument('-r',required=True, help='bed file of reference annotation')
    parser.add_argument('-o',required=True, help='bed file of isoforms')
    parser.add_argument('-tmp',type=str, default='/dev/shm/' , help='folder of tmp files')
    parser.add_argument('-bed',type=str, default='bedtools')
    parser.add_argument('-w',type= str,default="2.3", help='weight of intron')
    parser.add_argument('-t',type=str, default='1', help='number of threads')
    parser.add_argument('-s1',type=str, default='3', help="number of standard isoforms' support reads")
    parser.add_argument('-s2',type=str, default='5', help="number of high-confident isoforms' support reads")
    parser.add_argument('-q',required=True, help="file of read count of isoforms")
    args, unknown = parser.parse_known_args()

    intersect_command = [args.bed, 'intersect', '-wa', '-wb', '-r', '-a', args.i, '-b', args.r ]

    subprocess.call(intersect_command, stdout=open(args.tmp+'intersect_reads_ref.bed', 'w'))
    subprocess.call([sys.executable, 'bin/annotation_round1.py', args.tmp+'intersect_reads_ref.bed', args.i,args.r, args.tmp])
    subprocess.call([sys.executable, 'bin/round1_cluster.py', args.tmp + "/annotation_round1/", args.tmp + "/cluster_round1/",args.w, args.t, args.tmp + 'cluster/'])
    subprocess.call([sys.executable, 'bin/rep_round1.py',args.tmp, args.r, args.i])
    intersect_command = [args.bed, 'intersect', '-a',args.tmp+'novel.bed', '-b',args.tmp+'rep_round1.short.bed','-s', '-wao' ]
    subprocess.call(intersect_command, stdout=open(args.tmp+'/intersect_novel_rep_round1.bed', 'w'))
    sort_command = [args.bed, 'sort', '-i',args.tmp+'novel.bed']
    subprocess.call(sort_command, stdout=open(args.tmp+'novel.sort.bed','w'))
    merge_command = [args.bed,  'merge', '-c', '4', '-o', 'distinct', '-s', '-i', args.tmp+'novel.sort.bed']
    subprocess.call(merge_command, stdout=open(args.tmp+'merge_novelv.bed', 'w'))
    subprocess.call([sys.executable, 'bin/annotation_round2.py',args.tmp,args.i])
    subprocess.call([sys.executable, 'bin/round2_cluster.py',args.tmp+"/annotation_round2/", args.tmp+"/cluster_round2/",args.w,args.t,args.tmp+'cluster/'])
    subprocess.call([sys.executable, 'bin/rep_round2.py',args.tmp, args.s1])

    sort_command = [args.bed, 'sort', '-i',args.tmp+'rep_round2.short.bed']
    subprocess.call(sort_command, stdout=open(args.tmp+'rep_round2.sort.bed','w'))
    merge_command = [args.bed,  'merge', '-c', '4', '-o', 'distinct', '-s', '-i', args.tmp+'rep_round2.sort.bed']
    subprocess.call(merge_command, stdout=open(args.tmp+'merge_rep_round2.bed', 'w'))

    subprocess.call([sys.executable, 'bin/annotation_round3.py',args.tmp,args.i, args.r])
    subprocess.call([sys.executable, 'bin/round3_cluster.py',args.tmp+"/annotation_round3/", args.tmp+"/cluster_round3/",args.w,args.t,args.tmp+'cluster/'])
    subprocess.call([sys.executable, 'bin/support_reads.py',args.tmp,args.o,args.q,args.s2])
    rm_command = ['rm', '-r', args.tmp]
    subprocess.call(rm_command)
