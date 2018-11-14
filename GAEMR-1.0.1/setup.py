from distutils.core import setup
from os.path import join, dirname

setup(
    name='GAEMR',
    version='1.0.1',
    author='Genome Assembly and Analysis Group - Broad Institute',
    author_email='gaag@broadinstitute.org',
    packages=['gaemr'],
    scripts=[
	'bin/GAEMR.py',
	'bin/align_reads.py',
	'bin/analyze_gap_ends.py',
	'bin/analyze_rna_hits.py',
	'bin/assembly_stats_comparison.py',
	'bin/assess_readoids_coords.py',
	'bin/basic_assembly_stats.py',
	'bin/blast_bubbles.py',
	'bin/blast_map.py',
	'bin/compare_to_reference.py',
	'bin/create_readoids.py',
	'bin/gaemr_html.py',
	'bin/generate_bam_plots.py',
	'bin/get_bam_coverage_stats.py',
	'bin/get_blast_hit_taxonomy.py',
	'bin/get_simple_bam_stats.py',
	'bin/identify_coverage_anomalies.py',
	'bin/kmer_copy_number.py',
	'bin/kmer_coverage.py',
	'bin/make_detailed_assembly_table.py',
	'bin/make_standard_assembly_files.py',
	'bin/ono_to_ref.py',
	'bin/parse_blast_xml.py',
	'bin/plot_insert_size.py',
	'bin/read_format_converter.py',
	'bin/run_blast.py',
	'bin/run_contamination_screen.py',
	'bin/run_insert_size_from_bam.py',
	'bin/run_nucmer.py',
	'bin/run_rna_analysis.py',
	'bin/run_scaffold_accuracy.py',
	'bin/run_vecscreen.py',
	'bin/subset_fasta.py'
    ],
    url='http://www.broadinstitute.org/software/gaemr/',
    license='LICENSE.txt',
    description='Genome Assembly Evaluation Metrics and Reporting.',
    long_description=open(join(dirname(__file__), 'README.txt')).read(),
)