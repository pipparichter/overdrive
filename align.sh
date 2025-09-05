! mmseqs createdb ./overdrive.fn ./mmseqs/databases/overdriveDB --dbtype 2 

# Run prefilter with high sensitivity (-s 10) to matches. 
! mmseqs prefilter -s 7.5 ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB_prefilter

# ! mmseqs align --min-aln-len 8 ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB_prefilter ./mmseqs/overdriveDB_align
! mmseqs align ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB_prefilter ./mmseqs/overdriveDB_align

! mmseqs convertalis --search-type 2 --format-output "query,target,evalue,pident,bits,qseq,tseq,alnlen,qstart,qend,tstart,tend" ./mmseqs/databases/overdriveDB ./mmseqs/databases/overdriveDB ./mmseqs/overdriveDB_align ./overdrive_align.tsv
