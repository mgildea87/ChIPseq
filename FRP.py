import subprocess, pandas as pd

sample_file = "samples_info.tab"
sample = pd.read_table(sample_file)['Sample']
replicate = pd.read_table(sample_file)['Replicate']
condition = pd.read_table(sample_file)['Condition']

sample_ids = []
for i in range(len(sample)):
	sample_ids.append('%s_%s_%s' % (sample[i], condition[i], replicate[i]))
sample_ids = pd.unique(sample_ids).tolist()

total_fragments = []
peak_fragments = []
def FRP(sample):
	cmd_1 = ['bedtools', 'bamtobed', '-bedpe', '-i', 'alignment/%s_Antibody.bam' % (sample)] 
	cmd_2 = ['bedtools', 'intersect', '-a', 'stdin', '-b', 'peaks/%s_peaks.narrowPeak' % (sample)]
	cmd_3 = ['wc', '-l']
	step_1 = subprocess.Popen(cmd_1, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_2, stdin = step_1.stdout, stdout=subprocess.PIPE)
	step_3 = subprocess.Popen(cmd_3, stdin = step_2.stdout, stdout=subprocess.PIPE)
	frag_peaks = step_3.stdout.read()
	frag_peaks.strip()
	peak_fragments.append(float(frag_peaks))

	cmd_4 = ['samtools', 'view', 'alignment/%s_Antibody.bam' % (sample)]
	step_1 = subprocess.Popen(cmd_4, stdout=subprocess.PIPE)
	step_2 = subprocess.Popen(cmd_3, stdin = step_1.stdout, stdout=subprocess.PIPE)

	frag = step_2.stdout.read()
	frag.strip()
	frag = float(frag)
	total_fragments.append(frag/2)

for samples in sample_ids:
	FRP(samples)

with open('FRP.txt', 'w') as out:
	out.write('Sample\tTotal_fragments\tFragments_in_peaks\tFRP\n')
	for i in range(len(sample_ids)):
		out.write('%s\t%s\t%s\t%s\n' % (sample_ids[i], total_fragments[i], peak_fragments[i], peak_fragments[i]/total_fragments[i]))
	out.close()	
