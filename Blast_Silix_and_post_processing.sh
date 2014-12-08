# Pipeline for ALL-by-ALL BLAST and Silix clustering (e.g. for BovB):

# Step by step instructions for doing an ALL-by-ALL BLAST of hits, clustering them, then post-processing to look for potential HT candidates


# Make a FASTA file of all the hits (or 'cutoff' hits - the hits above the cutoff length threshold, that seem valid)
# Make sure every seq has a unique name and includes the species name
# Name this file 'BovB_cutoff_hits.fasta'

# Go on leeuwenhoek: /atma/TESTINGSILIX/

# Format the FASTA file of seqs as a database 
formatdb -i BovB_cutoff_hits.fasta -p F -o T -n BovBHits
# -p: set to F for nucleotide seqs
# -o: set to T to parse seqID and create indexes
# -n: basename for this formatted database, e.g. BovBHits

# Perform an all-against-all blast of the hits in the database
blastall -p blastn -d BovBHits -d BovBHits -i BovB_cutoff_hits.fasta -r 2 -F F -m 8 -e 1e-10 -a 4 -o blastall_BovB.out
# -p: chooses the BLAST program
# -d: database to query
# -i: query seq(s) - since we want to do an ALL-by-ALL comparison, the query seqs are the same as the database seqs
# -r: reward for nucleotide model. Setting this to 2 makes it more adapted for divergent seqs, facilitating the seq clustering by Silix
# -f: filter/mask query seq (set to False, since we're working with repeats)
# -m: for tabular output (to use with Silix), set to 8
# -e: e-value
# -a: number of CPUs (e.g. set to 4)
# -o: output file

# Move the output file blastall_BovB.out to computer: ~/atma/TESTINGSILIX/ (e.g. with Fugu)

# Use Silix software to cluster the hits at a range of different % ids
for j in 60 65 70 75 80 85 90 95 98
do
silix BovB_cutoff_hits.fasta blastall_BovB.out -f FAM -i "0.$j" -r 0.70 > seqBovB_$j.fnodes
done
# -f: prefix of families
# -i: min % identity
# -r: min length/overlap

# Sort the lines in each file numerically (e.g. FAM000001 to FAM-whatever-largest-num-is)
for j in 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.98
do
sort seqBovB_$j.fnodes > seqBovB_$j_sorted.fnodes
done