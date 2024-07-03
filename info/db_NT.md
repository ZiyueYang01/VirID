## NCBI Nucleotide database (NT) wihtout viruses
The following are the steps to exclude virus sequences from the NT database from the sequence file in the `.fasta` format ([NT fasta](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)).

### 1. Remove the viral sequences listed in [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/).
We have pre-collected all virus nucleic acid [sequence IDs](./NCBI_Virus_All_Nucleotides_Accession_List_20240617.7z) (as of 18 June 2024) available in the NCBI Virus website.


Example code
```shell
seqkit grep -f NCBI_Virus_DNA_viruses_AccessionList_20240617.tsv -v nt.fasta -o NT_noViruses_0.fas
```


### 2.Remove the viral sequences based on their taxonomic IDs (taxids)
[nucl_gb.accession2taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)

[taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)

```shell
seqkit seq -n -i NT_noViruses_0.fas > NT_noViruses_0.id

grep -F -f NT_noViruses_0.id nucl_gb.accession2taxid > NT_noViruses_0.id2taxid

sed -i '/^$/d' NT_noViruses_0.id2taxid


cut -f3 NT_noViruses_0.id2taxid \
| sort -u \
| taxonkit --data-dir ./taxdump lineage \
| awk '$2>0' \
| taxonkit --data-dir ./taxdump reformat -f "{k}\t{kID}\t{p}\t{pID}\t{c}\t{cID}\t{o}\t{oID}\t{f}\t{fID}\t{g}\t{gID}\t{s}\t{sID}" -F \
| cut -f1,2,4,6,8,10,12,14 \
> NT_noViruses_0.tax

awk -F'\t' '$2 != 10239' NT_noViruses_0.tax > NT_noViruses_1.tax

```
