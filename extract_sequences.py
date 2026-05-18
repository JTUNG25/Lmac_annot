#!/usr/bin/env python3
"""
Extract protein sequences from Leptosphaeria maculans JN3 annotation files
Supports GFF3/GTF formats — no gffutils or pandas required, only biopython + stdlib
"""

import os
import sys
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class GFFParser:
    """Lightweight GFF3 parser using only biopython/stdlib."""

    def __init__(self, gff_file):
        self.gff_file = gff_file
        # gene_id -> list of CDS dicts
        self.cds_by_gene = defaultdict(list)
        self._parse()

    def _get_attr(self, attr_string, key):
        """Extract a value from a GFF attribute string."""
        for part in attr_string.strip().split(";"):
            part = part.strip()
            if part.startswith(key + "="):
                return part[len(key) + 1 :]
            # GTF-style: key "value"
            if part.startswith(key + ' "'):
                return part[len(key) + 2 :].rstrip('"')
        return None

    def _parse(self):
        print(f"Parsing annotation file: {self.gff_file}")
        # Two-pass approach:
        #   Pass 1 — collect mRNA->gene mapping
        #   Pass 2 — collect CDS records, resolve to gene ID

        mrna_to_gene = {}

        with open(self.gff_file) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                ftype = cols[2]
                if ftype not in ("mRNA", "transcript"):
                    continue
                attrs = cols[8]
                mrna_id = self._get_attr(attrs, "ID")
                gene_id = self._get_attr(attrs, "Parent") or self._get_attr(
                    attrs, "gene_id"
                )
                if mrna_id and gene_id:
                    mrna_to_gene[mrna_id] = gene_id

        with open(self.gff_file) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                if cols[2] != "CDS":
                    continue

                seqid = cols[0]
                start = int(cols[3])  # 1-based inclusive
                end = int(cols[4])  # 1-based inclusive
                strand = cols[6]
                attrs = cols[8]

                parent = self._get_attr(attrs, "Parent") or self._get_attr(
                    attrs, "gene_id"
                )
                if not parent:
                    continue

                # Resolve to gene ID (parent may be mRNA or gene directly)
                gene_id = mrna_to_gene.get(parent, parent)

                self.cds_by_gene[gene_id].append(
                    {
                        "seqid": seqid,
                        "start": start,
                        "end": end,
                        "strand": strand,
                    }
                )

        print(f"Found CDS for {len(self.cds_by_gene)} genes")


class SequenceExtractor:
    def __init__(self, genome_fasta, annotation_file, output_dir="extracted_sequences"):
        self.genome_fasta = genome_fasta
        self.annotation_file = annotation_file
        self.output_dir = output_dir
        self.genome_dict = None
        os.makedirs(output_dir, exist_ok=True)

    def load_genome(self):
        print("Loading genome FASTA file...")
        try:
            self.genome_dict = SeqIO.to_dict(SeqIO.parse(self.genome_fasta, "fasta"))
            print(f"Loaded {len(self.genome_dict)} contigs/chromosomes")
        except Exception as e:
            print(f"Error loading genome: {e}")
            sys.exit(1)

    def _extract_cds_sequence(self, cds_list):
        """Concatenate CDS features into a single coding sequence."""
        strand = cds_list[0]["strand"]

        # Forward strand: ascending start; reverse strand: descending start
        sorted_cds = sorted(cds_list, key=lambda x: x["start"], reverse=(strand == "-"))

        pieces = []
        for cds in sorted_cds:
            contig = cds["seqid"]
            if contig not in self.genome_dict:
                print(f"  Warning: contig '{contig}' not in genome — skipping")
                continue
            # GFF is 1-based inclusive; Biopython slicing is 0-based half-open
            seq = self.genome_dict[contig].seq[cds["start"] - 1 : cds["end"]]
            if strand == "-":
                seq = seq.reverse_complement()
            pieces.append(seq)

        return sum(pieces, Seq("")) if pieces else None

    def extract_from_gff(self, gene_ids=None):
        """Extract protein and CDS records from GFF, optionally filtered to gene_ids."""
        parser = GFFParser(self.annotation_file)

        target_genes = set(gene_ids) if gene_ids else None

        genes_to_process = (
            {g: v for g, v in parser.cds_by_gene.items() if g in target_genes}
            if target_genes
            else parser.cds_by_gene
        )

        if target_genes:
            missing = sorted(target_genes - set(genes_to_process))
            if missing:
                print(f"Warning: {len(missing)} gene IDs not found in annotation:")
                for m in missing[:10]:
                    print(f"  {m}")
                if len(missing) > 10:
                    print(f"  ... and {len(missing) - 10} more")

        protein_records = []
        cds_records = []

        for gene_id, cds_list in genes_to_process.items():
            print(f"Processing gene: {gene_id}")
            cds_seq = self._extract_cds_sequence(cds_list)
            if not cds_seq:
                print(f"  Warning: no CDS sequence extracted for {gene_id}")
                continue

            cds_records.append(
                SeqRecord(
                    cds_seq, id=gene_id, description=f"CDS sequence for {gene_id}"
                )
            )

            try:
                prot = cds_seq.translate()
                if prot.endswith("*"):
                    prot = prot[:-1]
                protein_records.append(
                    SeqRecord(
                        prot, id=gene_id, description=f"Protein sequence for {gene_id}"
                    )
                )
            except Exception as e:
                print(f"  Warning: could not translate {gene_id}: {e}")

        print(
            f"\nExtracted {len(protein_records)} protein sequences "
            f"({len(cds_records)} CDS sequences)"
        )
        return protein_records, cds_records

    def extract_specific_genes(self, gene_list_file):
        print(f"Reading gene list from {gene_list_file}")
        with open(gene_list_file) as f:
            gene_ids = [line.strip() for line in f if line.strip()]
        print(f"Found {len(gene_ids)} gene IDs to extract")
        return self.extract_from_gff(gene_ids)

    def save_sequences(self, protein_records, cds_records):
        protein_file = os.path.join(self.output_dir, "proteins.fasta")
        cds_file = os.path.join(self.output_dir, "cds.fasta")

        if protein_records:
            SeqIO.write(protein_records, protein_file, "fasta")
            print(f"Saved {len(protein_records)} protein sequences → {protein_file}")
        else:
            open(protein_file, "w").close()  # empty file so Snakemake output exists

        if cds_records:
            SeqIO.write(cds_records, cds_file, "fasta")
            print(f"Saved {len(cds_records)} CDS sequences → {cds_file}")
        else:
            open(cds_file, "w").close()

        return protein_file, cds_file

    def create_gene_summary(self, protein_records, cds_records):
        """
        Write extraction_summary.csv using only stdlib csv module.

        csv.DictWriter works in two steps:
          1. writeheader() — writes the column names as the first row
          2. writerow({...}) — writes one data row per call, matching keys
             to column names. Called once per gene inside the for loop.
        """
        cds_len = {r.id: len(r.seq) for r in cds_records}
        summary_file = os.path.join(self.output_dir, "extraction_summary.csv")

        with open(summary_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["Gene_ID", "Protein_Length", "CDS_Length", "Description"]
            )
            writer.writeheader()  # row 1: Gene_ID,Protein_Length,CDS_Length,Description
            for r in protein_records:
                writer.writerow(
                    {  # one row per gene
                        "Gene_ID": r.id,
                        "Protein_Length": len(r.seq),
                        "CDS_Length": cds_len.get(r.id, 0),
                        "Description": r.description,
                    }
                )

        print(f"Saved gene summary → {summary_file}")
        return summary_file


def main():
    parser = argparse.ArgumentParser(
        description="Extract protein/CDS sequences from JN3 GFF3 annotation"
    )
    parser.add_argument("genome", help="Genome FASTA file")
    parser.add_argument("annotation", help="GFF3/GTF annotation file")
    parser.add_argument("--genes", help="File with gene IDs to extract (one per line)")
    parser.add_argument(
        "--output", default="extracted_sequences", help="Output directory"
    )
    args = parser.parse_args()

    extractor = SequenceExtractor(args.genome, args.annotation, args.output)
    extractor.load_genome()

    if args.genes:
        protein_records, cds_records = extractor.extract_specific_genes(args.genes)
    else:
        protein_records, cds_records = extractor.extract_from_gff()

    if not protein_records:
        print("No sequences extracted — check your genome, annotation, and gene IDs.")

    protein_file, cds_file = extractor.save_sequences(protein_records, cds_records)
    summary_file = extractor.create_gene_summary(protein_records, cds_records)

    print(f"\nExtraction complete!")
    print(f"  Proteins : {protein_file}")
    print(f"  CDS      : {cds_file}")
    print(f"  Summary  : {summary_file}")


if __name__ == "__main__":
    main()
