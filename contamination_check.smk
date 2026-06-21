#!/usr/bin/env python3

GENOME = "Lmac_D5"
ASSEMBLY = "data/genomes/Lmac_D5.fasta"
GFF3 = "results/tiberius_evidence/tiberius_train.gff3"
OMARK_RESULT = "results/omark/Lmac_D5_evidence/Lmac_D5_evidence.omamer"
OUTDIR = "results/contamination_screening"


rule target:
    input:
        f"{OUTDIR}/contamination_report.txt",
        f"{OUTDIR}/scaffold_summary.tsv",
        f"{OUTDIR}/species_breakdown.txt",


rule parse_omark:
    input:
        omark=OMARK_RESULT,
    output:
        contam_proteins=f"{OUTDIR}/contaminated_proteins_list.tsv",
        species_map=f"{OUTDIR}/protein_species_map.tsv",
    log:
        "logs/parse_omark.log",
    script:
        """
        from collections import defaultdict
        
        contaminated = {}
        all_proteins = {}
        target_species = {'Leptosphaeria maculans', 'Leptosphaeria', 'Lmac'}
        
        with open(input.omark) as f:
            for i, line in enumerate(f):
                if i == 0: continue
                fields = line.strip().split('\\t')
                if len(fields) < 3: continue
                
                query_id = fields[0]
                species_info = fields[2] if len(fields) > 2 else "UNKNOWN"
                all_proteins[query_id] = species_info
                
                is_target = any(target in species_info for target in target_species)
                if not is_target:
                    contaminated[query_id] = species_info
        
        with open(output.species_map, 'w') as f:
            f.write("protein_id\\tspecies\\tcontamination_status\\n")
            for protein_id in sorted(all_proteins.keys()):
                status = "CONTAMINATED" if protein_id in contaminated else "TARGET"
                f.write(f"{protein_id}\\t{all_proteins[protein_id]}\\t{status}\\n")
        
        with open(output.contam_proteins, 'w') as f:
            f.write("protein_id\\tspecies\\n")
            for protein_id in sorted(contaminated.keys()):
                f.write(f"{protein_id}\\t{contaminated[protein_id]}\\n")
        
        print(f"Total: {len(all_proteins)}, Contaminated: {len(contaminated)}")
        """


rule map_proteins_to_scaffolds:
    input:
        gff3=GFF3,
        assembly=ASSEMBLY,
        contam_proteins=f"{OUTDIR}/contaminated_proteins_list.tsv",
    output:
        scaffold_report=f"{OUTDIR}/contamination_report.txt",
        summary_tsv=f"{OUTDIR}/scaffold_summary.tsv",
    log:
        "logs/map_proteins_to_scaffolds.log",
    script:
        """
        from Bio import SeqIO
        from collections import defaultdict
        import pandas as pd
        
        # Load scaffold sizes
        scaffold_sizes = {}
        for record in SeqIO.parse(input.assembly, 'fasta'):
            scaffold_sizes[record.id] = len(record.seq)
        
        # Load contaminated proteins
        contaminated = set()
        with open(input.contam_proteins) as f:
            next(f)
            for line in f:
                contaminated.add(line.strip().split('\\t')[0])
        
        # Map proteins to scaffolds
        scaffold_proteins = defaultdict(list)
        with open(input.gff3) as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\\t')
                if len(fields) < 9 or fields[2] not in ['gene', 'mRNA']: continue
                
                scaffold = fields[0]
                attrs = fields[8]
                gene_id = None
                for attr in attrs.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('ID=')[1].split('.')[0]
                        break
                
                if gene_id:
                    is_contam = gene_id in contaminated
                    scaffold_proteins[scaffold].append((gene_id, is_contam))
        
        # Generate report
        scaffold_stats = []
        for scaffold in scaffold_proteins:
            proteins = scaffold_proteins[scaffold]
            total = len(proteins)
            num_contam = sum(1 for _, is_contam in proteins if is_contam)
            pct = 100 * num_contam / total if total > 0 else 0
            size = scaffold_sizes.get(scaffold, 0)
            
            scaffold_stats.append({
                'scaffold': scaffold,
                'size_bp': size,
                'total_proteins': total,
                'contaminated_proteins': num_contam,
                'pct_contaminated': pct,
            })
        
        scaffold_stats_sorted = sorted(scaffold_stats, key=lambda x: -x['pct_contaminated'])
        
        # Write report
        report_lines = ["CONTAMINATION ANALYSIS REPORT", "=" * 80, ""]
        report_lines.append(f"Total scaffolds: {len(scaffold_proteins)}")
        report_lines.append(f"Total contaminated proteins: {len(contaminated)}")
        report_lines.append("")
        
        high = [s for s in scaffold_stats_sorted if s['pct_contaminated'] >= 80]
        med = [s for s in scaffold_stats_sorted if 40 <= s['pct_contaminated'] < 80]
        low = [s for s in scaffold_stats_sorted if 0 < s['pct_contaminated'] < 40]
        
        report_lines.append(f"HIGH (≥80%): {len(high)}")
        report_lines.append(f"MEDIUM (40-80%): {len(med)}")
        report_lines.append(f"LOW (0-40%): {len(low)}")
        report_lines.append("")
        
        for s in scaffold_stats_sorted:
            if s['pct_contaminated'] > 0:
                report_lines.append(
                    f"{s['scaffold']:<40} {s['pct_contaminated']:>6.1f}% "
                    f"({s['contaminated_proteins']}/{s['total_proteins']} proteins)"
                )
        
        with open(output.scaffold_report, 'w') as f:
            f.write('\\n'.join(report_lines))
        
        # Write TSV
        df = pd.DataFrame(scaffold_stats_sorted)
        df = df.sort_values('pct_contaminated', ascending=False)
        df.to_csv(output.summary_tsv, sep='\\t', index=False)
        
        print(f"HIGH: {len(high)}, MED: {len(med)}, LOW: {len(low)}")
        """


rule species_breakdown:
    input:
        species_map=f"{OUTDIR}/protein_species_map.tsv",
    output:
        species_summary=f"{OUTDIR}/species_breakdown.txt",
    log:
        "logs/species_breakdown.log",
    script:
        """
        from collections import defaultdict
        
        species_count = defaultdict(int)
        status_count = defaultdict(int)
        
        with open(input.species_map) as f:
            next(f)
            for line in f:
                protein_id, species, status = line.strip().split('\\t')
                species_count[species] += 1
                status_count[status] += 1
        
        report_lines = ["SPECIES COMPOSITION SUMMARY", "=" * 80, ""]
        total = sum(status_count.values())
        
        for status, count in sorted(status_count.items(), key=lambda x: -x[1]):
            pct = 100 * count / total
            report_lines.append(f"  {status:<20} {count:>6} ({pct:>5.1f}%)")
        
        report_lines.append("")
        for species, count in sorted(species_count.items(), key=lambda x: -x[1])[:20]:
            pct = 100 * count / total
            is_target = 'TARGET' if 'Leptosphaeria' in species else 'CONTAMINANT'
            report_lines.append(f"  {species:<50} {count:>6} ({pct:>5.1f}%) [{is_target}]")
        
        with open(output.species_summary, 'w') as f:
            f.write('\\n'.join(report_lines))
        
        print('\\n'.join(report_lines))
        """
