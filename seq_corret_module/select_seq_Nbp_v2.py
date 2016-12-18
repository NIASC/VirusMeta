import os, sys
from Bio import SeqIO

def select_seq_Nbp(min_bp,max_bp,in_handle,out_handle):
    '''
    select sequences according to specified read length
    and rename them according to file name
    '''
    query_fasta_name = os.path.split(os.path.abspath(in_handle))[1] #query fasta name
    #query_fasta_name = query_fasta_name.replace (".","_")
    query_fasta_name = query_fasta_name.replace (".","@")
    query_fasta_name = query_fasta_name.split ("@")[0]
    print query_fasta_name
    records  = SeqIO.parse(open(in_handle, 'rU'), 'fasta')
    if os.path.exists(out_handle):
       out_handle = open(out_handle, 'a')
    else:
       out_handle = open(out_handle, 'w')
    seq_nr = 0   
    for record in records:
        if len(record.seq) >= int(min_bp) and len(record.seq) <= int(max_bp):
           #print len(record.seq)
           out_handle.write(">seq%i_%s\n%s\n" % (seq_nr,query_fasta_name,record.seq))
           seq_nr += 1


#select_seq_Nbp(200,1000,"%s/Trinity.fasta")
select_seq_Nbp(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

