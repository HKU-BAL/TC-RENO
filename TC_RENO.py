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
    parser.add_argument('-s1',type=str, default='6', help="number of isoforms' support reads")
    parser.add_argument('-s2',type=str, default='10', help="number of short isoforms' support reads")
    parser.add_argument('-q',required=True, help="file of read count of isoforms")

    args = parser.parse_args() 
    args, unknown = parser.parse_known_args()


    intersect_command = [args.bed, 'intersect', '-s', '-wao', '-e', '-F', '0.9','-f','0.9', '-a', args.i, '-b', args.r ]
    subprocess.call(intersect_command, stdout=open(args.tmp+'/intersect_reads_gene.bed', 'w'))
    subprocess.call([sys.executable, 'bin/annotation_round1.py', args.tmp,args.i])
    subprocess.call([sys.executable, 'bin/cluster1.py',args.tmp+"/annotation_round1/", args.tmp+"/cluster_round1_TC_RENO/",args.w,args.t,args.tmp+'/cluster_tmp/'])
    subprocess.call([sys.executable, 'bin/rep_round1.py', args.i, args.tmp])
    sort_command = [args.bed, 'sort', '-i',args.tmp+'/rep_round1.bed']
    subprocess.call(sort_command, stdout=open(args.tmp+'/rep_round1.sort.bed','w'))
    merge_command = [args.bed,  'merge', '-c', '4', '-o', 'distinct', '-s', '-i', args.tmp+'/rep_round1.sort.bed']
    subprocess.call(merge_command, stdout=open(args.tmp+'/merge_round1.bed', 'w'))
    subprocess.call([sys.executable, 'bin/annotation_round2.py',args.tmp,args.i])
    subprocess.call([sys.executable, 'bin/cluster1.py',args.tmp+"/annotation_round2/", args.tmp+"/cluster_round2_TC_RENO/",args.w,args.t,args.tmp+'/cluster_tmp/'])
    subprocess.call([sys.executable, 'bin/rep_round2.py',args.tmp,args.s1])
    subprocess.call([sys.executable, 'bin/cluster2.py',args.tmp+"/annotation_round3/", args.tmp+"/cluster_round3_TC_RENO/",args.w,args.t,args.tmp+'/cluster_tmp/'])
    subprocess.call([sys.executable, 'bin/support_reads.py',args.tmp,args.o,args.q,args.s1, args.s2, args.i])

    rm_command = ['rm', '-r', args.tmp]
    subprocess.call(rm_command)
