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
silix BovB_cutoff_hits.fasta blastall_BovB.out -i "0.$j" -r 0.70 > seqBovB_$j.fnodes
done
# -f: prefix of families - don't need to set a prefix!
# -i: min % identity
# -r: min length/overlap

# Sort the lines in each file numerically (e.g. family 1 to whatever-largest-num-is)
for j in 60 65 70 75 80 85 90 95 98
do
sort -n seqBovB_$j.fnodes > seqBovB_sorted_$j.fnodes # Can't name output file seqBovB_$j_sorted.fnodes or anything with $j_*.fnodes - causes overwriting error 
done

# Separate each family/cluster by a blank line
for j in 60 65 70 75 80 85 90 95 98
do
awk -v i=1 'NR>1 && $i!=p { print "" }{ p = $i } 1' seqBovB_sorted_$j.fnodes > seqBovB_separated_$j.fnodes
done
# i.e. on any line after the first, if the value of the "i"-th column (where i=1) is different to the previous value, print a blank line. 
# Always set the value of p
# 1 at the end means true - so that awk prints the line

# Remove singleton clusters
# Haven't figured out how to do this with 'awk' yet
# For now, open the files in TextWrangler and delete the singletons (at the bottom, after the clusters)

# Similarly (also need to figure out), remove species-specific (monophyletic) clusters, to leave only clusters containing >1 species

# Remove clusters than seem to be composed of (really) short fragments only (e.g. 200bp)

# Want to separate 2nd column (seq) into multiple columns by replacing all '_' with tabs
# Most seqs are of the form Animal_seq_startCoord_endCoord (+ maybe.rev at the end)
# _seq_ may be chr, gi533376418embHF963279, anything so long as it is continous (i.e. not separated by '_')
# This means that, once separated, the 3rd field will be the startCoord
# But some genomes don't match this form - Silkworm, RockHyrax, TasDevil, Megabat, Elephant
# (e.g. Silkworm_gi509958623refNW_004586076.1_18_428; RockHyrax_scaffold_6197_26521_29396; TasDevil_chrX_GL867900_random_129412_131616; Megabat_scaffold_9013_11334_14329.rev; Elephant_scaffold_9_45655231_45658403)
# Need to alter these genomes to fit the required form (these must be reversible changes!!)
# Silkworm - Find: refNW_, Replace: refNW
# RockHyrax - Find: RockHyrax_scaffold_, Replace: RockHyrax_scaffold
# TasDevil - Find: _random_, Replace: random_
# Megabat - Find: Megabat_scaffold_, Replace: Megabat_scaffold
# Elephant - Find: Elephant_scaffold_, Replace: Elephant_scaffold
# Then replace all '_' with tabs - Find: _, Replace: tab (note: need to copy and paste a tab for it to work)
# Then for those files with '.rev' attached - Find: .rev, Replace: \t.rev
# Save changes as new file: e.g. seqBovB_tabdelim_60.fnodes

# Need to subtract 3000 from the 4th column, and add 3000 to the 5th column (for Extended Hits)
awk '{ print $1 "\t" $2 "\t" $3 "\t" ($4 - 3000) "\t" ($5 + 3000) "\t" $6 }' seqBovB_tabdelim_60.fnodes > seqBovB_extended_60.fnodes



