#!/bin/bash
"""
Complete setup script for L. maculans annotation on Bunya HPC
This script creates all necessary files and directories
"""

echo "=================================================="
echo "Setting up L. maculans Annotation Pipeline on HPC"
echo "Location: $(pwd)"
echo "=================================================="

# Create directory structure
mkdir -p {data/genome,data/gene_lists,scripts,logs,results,databases,profiles/bunya}

# Create the Bunya profile for Snakemake
cat > profiles/bunya/config.yaml << 'EOF'
cluster: 
  mkdir -p logs/{rule} &&
  sbatch
    --account=a_qaafi_chs
    --partition=general
    --qos=normal
    --job-name=smk-{rule}-{wildcards}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}M
    --time={resources.runtime}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
    --error=logs/{rule}/{rule}-{wildcards}-%j.err
default-resources:
  - mem_mb=8000
  - runtime=120
restart-times: 2
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
use-singularity: true
singularity-args: "--cleanenv --containall"
EOF

# Create the main Snakefile
cat > annotation.smk << 'EOF'
#!/usr/bin/env python3

# Container definitions
eggnog = "docker://quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
interproscan = "docker://quay.io/biocontainers/interproscan:5.62_94.0--hec16e2b_1"
diamond = "docker://quay.io/biocontainers/diamond:2.1.8--h43eeafb_0"
biopython = "docker://quay.io/biocontainers/biopython:1.81"

# Input files - EDIT THESE PATHS
GENOME = "data/genome/JN3.fasta"
ANNOTATION = "data/genome/JN3_annotation.gff3"

# Mutant gene lists - EDIT THESE
MUTANTS = {
    "rnai_targets": "data/gene_lists/rnai_genes.txt",
}

# Pipeline settings
RUN_INTERPROSCAN = False  # Set to True for full analysis (adds 3+ hours)
MAX_EVALUE = 1e-5

# Get all mutant names
mutant_names = list(MUTANTS.keys())

rule target:
    input:
        expand("results/annotation/{mutant}/final_annotation.xlsx", mutant=mutant_names),
        expand("results/annotation/{mutant}/readable_summary.csv", mutant=mutant_names)

rule extract_sequences:
    input:
        genome=GENOME,
        annotation=ANNOTATION,
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        proteins="results/sequences/{mutant}/proteins.fasta",
        cds="results/sequences/{mutant}/cds.fasta",
        summary="results/sequences/{mutant}/extraction_summary.csv"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30,
    container:
        biopython
    log:
        "logs/extract_sequences/{mutant}.log"
    shell:
        """
        mkdir -p results/sequences/{wildcards.mutant}
        
        python scripts/extract_sequences.py {input.genome} {input.annotation} \
            --genes {input.gene_list} \
            --output results/sequences/{wildcards.mutant} \
            > {log} 2>&1
        """

rule eggnog_annotation:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta"
    output:
        annotations="results/eggnog/{mutant}/eggnog.emapper.annotations",
        hits="results/eggnog/{mutant}/eggnog.emapper.hits",
        seed_orthologs="results/eggnog/{mutant}/eggnog.emapper.seed_orthologs"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120,
    container:
        eggnog
    log:
        "logs/eggnog/{mutant}.log"
    shell:
        """
        mkdir -p results/eggnog/{wildcards.mutant}
        
        emapper.py \
            -i {input.proteins} \
            --output results/eggnog/{wildcards.mutant}/eggnog \
            --output_dir results/eggnog/{wildcards.mutant} \
            -m diamond \
            --tax_scope fungi \
            --go_evidence non-electronic \
            --target_orthologs all \
            --seed_ortholog_evalue 0.001 \
            --seed_ortholog_score 60 \
            --override \
            --cpu {threads} \
            > {log} 2>&1
        """

rule diamond_blast:
    input:
        proteins="results/sequences/{mutant}/proteins.fasta",
        database="databases/uniprot_sprot.dmnd"
    output:
        results="results/diamond/{mutant}/blast_results.tsv"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
    container:
        diamond
    log:
        "logs/diamond/{mutant}.log"
    shell:
        """
        mkdir -p results/diamond/{wildcards.mutant}
        
        diamond blastp \
            -d {input.database} \
            -q {input.proteins} \
            -o {output.results} \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
            --max-target-seqs 5 \
            --evalue {MAX_EVALUE} \
            --threads {threads} \
            > {log} 2>&1
        """

rule combine_annotations:
    input:
        eggnog="results/eggnog/{mutant}/eggnog.emapper.annotations",
        diamond="results/diamond/{mutant}/blast_results.tsv",
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        summary="results/annotation/{mutant}/annotation_summary.csv"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    container:
        biopython
    log:
        "logs/combine_annotations/{mutant}.log"
    script:
        "scripts/combine_annotations.py"

rule final_annotation:
    input:
        combined="results/annotation/{mutant}/combined_annotations.xlsx",
        gene_list=lambda wildcards: MUTANTS[wildcards.mutant]
    output:
        final="results/annotation/{mutant}/final_annotation.xlsx",
        readable="results/annotation/{mutant}/readable_summary.csv"
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    container:
        biopython
    log:
        "logs/final_annotation/{mutant}.log"
    script:
        "scripts/create_final_report.py"

rule setup_databases:
    output:
        uniprot_db="databases/uniprot_sprot.dmnd"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    container:
        diamond
    log:
        "logs/setup_databases.log"
    shell:
        """
        mkdir -p databases
        cd databases
        
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip uniprot_sprot.fasta.gz
        
        diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
        rm uniprot_sprot.fasta
        
        echo "Database setup complete" > ../logs/setup_databases.log
        """
EOF

# Create simplified Python scripts for Snakemake
cat > scripts/extract_sequences.py << 'EOF'
#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO
import gffutils
import argparse
from pathlib import Path

def extract_sequences(genome_file, gff_file, gene_list_file, output_dir):
    """Extract protein sequences for specific genes"""
    
    # Read gene list
    with open(gene_list_file, 'r') as f:
        target_genes = set(line.strip() for line in f if line.strip())
    
    print(f"Extracting sequences for {len(target_genes)} genes")
    
    # Load genome
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    # Create GFF database
    db = gffutils.create_db(gff_file, ':memory:', force=True, keep_order=True)
    
    protein_records = []
    cds_records = []
    
    for gene_id in target_genes:
        try:
            gene = db[gene_id]
            
            # Find CDS features
            cds_features = list(db.children(gene, featuretype='CDS'))
            
            if cds_features:
                # Extract CDS sequence
                cds_seq = ""
                for cds in sorted(cds_features, key=lambda x: x.start):
                    contig_id = cds.seqid
                    if contig_id in genome_dict:
                        seq = genome_dict[contig_id].seq[cds.start-1:cds.end]
                        if cds.strand == '-':
                            seq = seq.reverse_complement()
                        cds_seq += str(seq)
                
                if cds_seq:
                    # Create CDS record
                    cds_record = SeqIO.SeqRecord(
                        seq=cds_seq,
                        id=gene_id,
                        description=f"CDS for {gene_id}"
                    )
                    cds_records.append(cds_record)
                    
                    # Translate to protein
                    try:
                        protein_seq = cds_record.seq.translate()
                        if protein_seq.endswith('*'):
                            protein_seq = protein_seq[:-1]
                        
                        protein_record = SeqIO.SeqRecord(
                            seq=protein_seq,
                            id=gene_id,
                            description=f"Protein for {gene_id}"
                        )
                        protein_records.append(protein_record)
                    except:
                        print(f"Warning: Could not translate {gene_id}")
        except:
            print(f"Warning: Could not find {gene_id}")
    
    # Save sequences
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    protein_file = output_dir / "proteins.fasta"
    cds_file = output_dir / "cds.fasta"
    
    SeqIO.write(protein_records, protein_file, "fasta")
    SeqIO.write(cds_records, cds_file, "fasta")
    
    print(f"Saved {len(protein_records)} protein sequences to {protein_file}")
    print(f"Saved {len(cds_records)} CDS sequences to {cds_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genome", help="Genome FASTA file")
    parser.add_argument("annotation", help="GFF3 annotation file")
    parser.add_argument("--genes", required=True, help="Gene list file")
    parser.add_argument("--output", required=True, help="Output directory")
    
    args = parser.parse_args()
    extract_sequences(args.genome, args.annotation, args.genes, args.output)
EOF

# Create other required Python scripts
cat > scripts/combine_annotations.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd
import sys
from pathlib import Path

# Simple annotation combination for Snakemake
eggnog_file = snakemake.input.eggnog
diamond_file = snakemake.input.diamond
output_combined = snakemake.output.combined

# Read eggNOG results
try:
    eggnog_df = pd.read_csv(eggnog_file, sep='\t', comment='#')
    eggnog_df = eggnog_df[['query', 'evalue', 'score', 'Description', 'Preferred_name', 'GOs', 'KEGG_Pathway']].copy()
    eggnog_df.rename(columns={'query': 'Gene_ID', 'Description': 'Functional_Description'}, inplace=True)
except:
    eggnog_df = pd.DataFrame()

# Read DIAMOND results
try:
    diamond_df = pd.read_csv(diamond_file, sep='\t', header=None, 
                           names=['Gene_ID', 'Subject_ID', 'Identity', 'Length', 'Mismatches',
                                  'Gap_Opens', 'Query_Start', 'Query_End', 'Subject_Start', 
                                  'Subject_End', 'E_value', 'Bit_Score', 'Subject_Description'])
    diamond_df = diamond_df.loc[diamond_df.groupby('Gene_ID')['E_value'].idxmin()]
    diamond_df = diamond_df[['Gene_ID', 'Identity', 'E_value', 'Subject_Description']].copy()
except:
    diamond_df = pd.DataFrame()

# Combine annotations
if not eggnog_df.empty and not diamond_df.empty:
    combined_df = eggnog_df.merge(diamond_df, on='Gene_ID', how='outer')
elif not eggnog_df.empty:
    combined_df = eggnog_df
elif not diamond_df.empty:
    combined_df = diamond_df
else:
    combined_df = pd.DataFrame()

# Save results
combined_df.to_excel(output_combined, index=False)

# Summary
summary_stats = {'Total_Genes': len(combined_df)}
pd.DataFrame([summary_stats]).to_csv(snakemake.output.summary, index=False)
EOF

cat > scripts/create_final_report.py << 'EOF'
#!/usr/bin/env python3
import pandas as pd

# Read combined annotations
df = pd.read_excel(snakemake.input.combined)

# Create readable summary
readable_df = df.copy()
readable_df['Confidence'] = 'Medium'  # Simplified confidence

# Save results
readable_df.to_excel(snakemake.output.final, index=False)
readable_df.to_csv(snakemake.output.readable, index=False)
EOF

# Create submission script
cat > submit_annotation.sh << 'EOF'
#!/bin/bash
#SBATCH --account=a_qaafi_chs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=lmac_annotation
#SBATCH --output=annot.log

source /sw/local/rocky8/noarch/rcc/software/miniforge/24.11.3-0/etc/profile.d/conda.sh
conda activate snakemake8
cd /QRISdata/Q9140/lmac/annot_lmac
export TMPDIR=$HOME/tmp
mkdir -p $TMPDIR

echo "=================================================="
echo "L. maculans Target-Specific Annotation Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=================================================="

# Check if databases exist, setup if needed
if [ ! -f "databases/uniprot_sprot.dmnd" ]; then
    echo "Setting up UniProt database (first time only)..."
    snakemake -s annotation.smk --profile profiles/bunya/ setup_databases
fi

if [ ! -d "databases/eggnog" ]; then
    echo "Setting up eggNOG database (first time only)..."
    conda activate eggnog-mapper
    download_eggnog_data.py --data_dir databases/eggnog
    conda activate snakemake8
fi

echo "✓ All databases ready"

# Run annotation pipeline
echo "Running annotation pipeline..."
snakemake -s annotation.smk --profile profiles/bunya/

echo "=================================================="
echo "Pipeline completed!"
echo "End time: $(date)"
echo "=================================================="

# Print results
echo "=== RESULTS ==="
find results -name "*.xlsx" -o -name "*.csv" | head -10
EOF

# Create example gene list
cat > data/gene_lists/rnai_genes.txt << 'EOF'
Lmb_jn3_04618
Lmb_jn3_04619
Lmb_jn3_04620
Lmb_jn3_04621
Lmb_jn3_04622
EOF

# Make scripts executable
chmod +x submit_annotation.sh
chmod +x scripts/*.py

echo ""
echo "=================================================="
echo "HPC Setup Complete!"
echo "=================================================="
echo ""
echo "Next steps:"
echo "1. Create eggnog-mapper environment:"
echo "   conda create -n eggnog-mapper -c bioconda eggnog-mapper"
echo ""
echo "2. Copy your data files:"
echo "   cp /path/to/JN3.fasta data/genome/"
echo "   cp /path/to/JN3_annotation.gff3 data/genome/"
echo ""
echo "3. Edit your gene list:"
echo "   nano data/gene_lists/rnai_genes.txt"
echo ""
echo "4. Submit job:"
echo "   sbatch submit_annotation.sh"
echo ""
echo "5. Monitor progress:"
echo "   squeue -u \$USER"
echo "   tail -f annot.log"
echo ""