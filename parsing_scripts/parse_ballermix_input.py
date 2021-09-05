'''
This script jointly reads the recombination map (hapmap format) and alignment file (.axt or .axt.gz) to:
	- polarize allele frequency
	- match recombination frequency
and - write outputs that matches the input format required by BalLeRMix

Note that this script has only been tested on human data. Please exercise caution when applied to other datasets, especially if the organism is not diploid. 
'''
import sys, os, re, gzip

SNP = {'A','T','C','G'}
int_regex = re.compile(r'[0-9]+')

def _get_ids_from_vcf_header(header, pop_list):
	# prep header
	header = header.strip().split('\t')
	#0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
	if pop_list is not None:
		# read pop_list file and remember the IDs
		ID_file = open(pop_list,'r')
		ID_list = ID_file.read().strip().split(",")
		ID_file.close()
		# get ID list
		ID_list = set(map( header.index , ID_list ))
		# sanity check
		assert min(ID_list) >= 9
	else:
		ID_list = range(9, len(header))

	print(f'Data from {len(ID_list)} samples will be counted. The total number of samples in this vcf is {len(header[9:])}.')

	return(ID_list)

def _get_GT_index(line):
	# get GT index
	GT_index = line[8].split(":").index("GT")
	return(GT_index)


def _open_large_file(filename, suffix = ""):
	if filename.lower().endswith('.gz'):
		file =  gzip.open(filename,'rt')
	elif filename.lower().endswith(suffix):
		file = open(filename, 'r')
	else:
		print(f'Unrecognized {suffix.upper()} file name. Please make sure it\'s in {suffix} or {suffix}.gz format.')
		sys.exit(1)

	return(file)


'''In the absence of rec map and alignment, read vcf and write MAF output with uniform recombination rate.'''
def parse_vcf_only(chrom, vcffile, rec_rate, outfile, pop_list):
	# write header of output
	out = open(outfile, 'w')
	out.write(f'position\tgenPos\tx\tn\n')

	# open vcf
	vcf = _open_large_file(vcffile, '.vcf')

	# start reading vcf
	l = vcf.readline()
	## go to the header
	while l.startswith('##'):
		l = vcf.readline()
		continue

	## now l is vcf header
	## from "FORMAT" get how to obtain GT
	individual_ids = _get_ids_from_vcf_header(l, pop_list)
	# start reading data
	#int_regex = re.compile(r'[0-9]+')
	#count_regex = re.compile(r'(0|1)')
	#0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
	l = vcf.readline()
	while l != '':
		l = l.strip().split('\t')
		if l[0] not in {chrom, 'chr'+chrom}:
			continue
		# only process bi-allelic SNPs
		if l[3] in SNP and l[4] in SNP and "PASS" in l[6]:
			pos = l[1]
			GT_id = _get_GT_index(l)
			x = 0; n = 0
			for i in individual_ids:
				GT = l[i].split(":")[GT_id]
				count = tuple(map( int, int_regex.findall(GT) ))
				# sanity check
				try:
					assert max(count) <= 1
				except Exception as e:
					print("Warning: This script only applies to diploids.")
					print(e)
					sys.exit(1)
				x += sum(count)
				n += len(count)
			# write output
			if x > 0 and x < n:
				x = min(n-x, x) #make sure it's MAF
				out.write(f"{pos}\t{float(pos) * rec_rate}\t{x}\t{n}\n")
		# read next line
		l = vcf.readline()

	# end of file. Close both input and output
	vcf.close()
	out.close()


'''When recombination map is given, jointly read vcf and bedfile to write MAF output'''
def parse_vcf_w_recmap(chrom, vcffile, rec_map, rec_rate, outfile, pop_list):
	# write header of output
	out = open(outfile, 'w')
	out.write(f'position\tgenPos\tx\tn\n')

	# open vcf
	vcf = _open_large_file(vcffile, '.vcf')

	# open recomb map
	recmap = _open_large_file(rec_map)
	## start reading map:
	rec_l =  recmap.readline().strip().split("\t")
	# skip header if there's any
	if int_regex.findall(rec_l[1]) != rec_l[1]:
		print('Skipping header.')
		rec_l =  recmap.readline().strip().split("\t")
	rec_ch = rec_l[0]
	assert rec_ch in {chrom, 'chr'+chrom }
	r_pos, rate, cum_cM = map(float, rec_l[1:])
	if cum_cM == 0:
		rate = 0
	elif rate == 'NA':
		rate = rec_rate * 1e6 #typically in unit of cM/Mb
	last_r_pos, last_r_rate, last_cum_CM = 0, rate, 0
	last_pos, last_gen_pos = 0, 0

	# start reading vcf
	l = vcf.readline()
	## go to the header
	while l.startswith('##'):
		l = vcf.readline()
		continue
	## now l is vcf header
	## from "FORMAT" get how to obtain GT
	individual_ids = _get_ids_from_vcf_header(l, pop_list)

	# start reading data
	#0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
	l = vcf.readline()	
	last_pos, last_gen_pos = 0, 0
	while l != '':
		l = l.strip().split('\t')
		# sanity check: make sure vcf and bed are on the same chromosome
		try:
			assert l[0] in {chrom, 'chr'+chrom }
			assert int_regex.findall(chrom)[0] == int_regex.findall(l[0])[0]
		except Exception as e:
			#print(ch, l[:5], type(ch), type(l[0]))
			print('Please make sure the recombination map and vcf cover the same chromosome.')
			print(l[0], e)
			sys.exit(1)

		# only process bi-allelic SNPs
		if l[3] in SNP and l[4] in SNP and "PASS" in l[6]:
			pos = int(l[1])
			# get genetic position (in cM)
			if pos <= r_pos:
				# use the rec rate in last seg in recmap
				gen_pos = last_gen_pos + (pos - last_pos) * (last_r_rate)/1e6
				# update
				last_pos, last_gen_pos = pos, gen_pos
			# otherwise, read another line
			else:
				while pos > r_pos:
					# update and read another line from recmap
					last_r_pos, last_r_rate, last_cum_CM = r_pos, rate, cum_cM
					rec_l =  recmap.readline()
					rec_l = rec_l.strip().split("\t")
					if rec_l == ['']: #presumably end-of-file
						print(f'End-of-file for the recombination map. Assume uniform rate of {rec_rate} after {last_r_pos}.')
						r_pos, rate, cum_cM = 8e8, rec_rate*1e6 , (last_cum_CM + rec_rate*(8e8-last_r_pos))
						break
					rec_ch = rec_l[0]
					try:
						assert rec_ch in {chrom , 'chr'+chrom}
					except Exception as e:
						print(e)
						print(rec_ch, chrom, rec_l)
					r_pos, rate, cum_cM = map(float, rec_l[1:])
					if rate == 'NA':
						rate = rec_rate * 1e6 #typically in unit of cM/Mb
				# now pos <= r_pos
				if last_pos < last_r_pos:
					gen_pos = last_cum_CM + (pos - last_r_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
				else:
					gen_pos = last_gen_pos + (pos - last_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
				# update gen_pos
				last_pos, last_gen_pos = pos, gen_pos

			x = 0; n = 0
			GT_id = _get_GT_index(l)
			for i in individual_ids:
				GT = l[i].split(":")[GT_id]
				count = tuple(map( int, int_regex.findall(GT) ))
				# sanity check
				try:
					assert max(count) <= 1
				except Exception as e:
					print("Warning: This script only applies to diploids.")
					print(e)
					sys.exit(1)
				x += sum(count)
				n += len(count)
			# write output if x != 0
			if x > 0 and x < n:
				x = min(n-x, x)
				out.write(f"{pos}\t{gen_pos}\t{x}\t{n}\n")
		# read next line
		l = vcf.readline()
	# end of file. Close all input and output
	vcf.close()
	out.close()
	recmap.close()

'''Read & "remember" the ancestral sequence on a chromosome
# make sure that chrom is chromosome id (without chr)'''
def _load_alignment(axtfile, chrom):
	axt_lib = {} #position: (ref, anc)
	aligned_positions = []
	atcg_start = re.compile(r'^[A|T|C|G]')
	#lone_atcg = re.compile(r'^[A|T|C|G]$')
	good_chr =  re.compile(r'^((chr)|(ch))?([0-9]+|X|Y)[a-z|A-Z]*$')
	axt = _open_large_file(axtfile, '.axt')
	l = axt.readline() 
	#skip annotations
	while l.startswith('#'):
		l = axt.readline()
	# first summary line
	assert not atcg_start.match(l)
	#print(l)
	#format: 0# 1chr 2start 3end 4chr_anc 5start_anc 6end_anc 7strand 8blastscore
	chroms = []
	while l != '':
		l = l.strip().split(" ")
		#print(l)
		#skip the next 4 lines if it's not for the current chr or the matching sequence is not from the legit assembly
		while l[1] not in {chrom, 'chr'+chrom} or not good_chr.match(l[4]):#
			chroms.append((l[1], l[4]));
			l = axt.readline() 
			l = axt.readline() 
			l = axt.readline() 
			l = axt.readline().strip().split(" ")
			if l == ['']:
				break
		if l == ['']:
			#print('skipped chromosomes:', set(chroms), '\n', chrom, 'chr'+chrom )
			break
		#l = l.strip().split(" ")
		start, end = int(l[2]), int(l[3])
		primary_seq = axt.readline().strip().upper()
		aligning_seq = axt.readline().strip().upper()
		if l == '':
			break
		pos = start - 1
		for i in range(len(primary_seq)):
			#if lone_atcg.match(primary_seq[i]):
			if primary_seq[i] in SNP:
				pos += 1
				# only record when there's no gap
				#if lone_atcg.match(aligning_seq[i]):
				if aligning_seq[i] in SNP:
					#print(f'recording {pos}')
					aligned_positions.append(pos)
					axt_lib[pos] = (primary_seq[i], aligning_seq[i])
		#read the next alignment
		l = axt.readline() #"\n"
		l = axt.readline() # next summary line
	#close file
	axt.close()
	#make sure list of positions are sorted
	aligned_positions.sort()
	# out of curiosity
	print(f'axt_lib takes memory {sys.getsizeof(axt_lib)}; aligned_pos list takes memory of {sys.getsizeof(aligned_positions)}.')
	#return lib
	return(axt_lib, aligned_positions)


'''Read through the sequence alignment, polarize the frequency, and get output, assume uniform rec rate'''
def parse_axt_and_vcf(ch, axtfile, vcffile, rec_rate, outfile, pop_list):
	# load axt file
	print(f'Loading the alignment for chr{ch}...')
	axt_lib, aligned_pos = _load_alignment(axtfile, ch)

	# read and parse vcf:
	# write header of output
	out = open(outfile, 'w')
	out.write(f'position\tgenPos\tx\tn\n')

	# open vcf
	vcf = _open_large_file(vcffile, '.vcf')

	# start reading vcf
	l = vcf.readline()
	## go to the header
	while l.startswith('##'):
		l = vcf.readline()
		continue
	## now l is vcf header
	## from "FORMAT" get how to obtain GT
	individual_ids = _get_ids_from_vcf_header(l, pop_list)
	try:
		sampleSize = 2* len(individual_ids)
	except Exception as e:
		print(e)
		print(individual_ids)
	# start reading data
	#int_regex = re.compile(r'[0-9]+')
	#0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
	l = vcf.readline()
	l = vcf.readline()
	try:
	#initiate indexing on aligned pos list
		axt_pos_i = 0 ; axt_pos = aligned_pos[axt_pos_i] ; axt_eof = False
	except Exception as e:
		print(e)
		print(aligned_pos, len(axt_lib))

	# start reading data
	while l != ['']:
		l = vcf.readline()
		if l == '':
			break
		l = l.strip().split('\t')
		# sanity check: it's the same chromosome
		assert l[0] in {ch, 'chr'+ch} 
		# check if pos in axt_pos & update axt_pos_i
		pos = int(l[1])
		while axt_pos < pos:
			axt_ref, axt_anc = axt_lib[axt_pos]
			# assume axt_ref; record as x=n=sampleSize
			if axt_ref != axt_anc:
				# write output
				out.write(f"{axt_pos}\t{axt_pos * rec_rate}\t{sampleSize}\t{sampleSize}\n")
			axt_pos_i += 1
			try:
				axt_pos = aligned_pos[axt_pos_i]
			except IndexError:
				print(f'No more recorded positions.')
				axt_eof = True
				break
		# end of while-loop for axt
		if axt_eof:
			break
		## Now axt_pos >= pos
		# only process bi-allelic SNPs
		if l[3] in SNP and l[4] in SNP and "PASS" in l[6]:
			# skip if not mapped
			if pos not in axt_lib:
				continue
			# skip if not bi-allelic
			elif len(set([axt_lib[pos], l[3], l[4]] )) > 2:
				continue
			# polarize
			else:
				assert pos == axt_pos
				axt_ref, axt_anc = axt_lib[pos]
				ref, alt = l[3], l[4]
				# if the alternate allele is ancestral, current count is for derived allele
				if ref == axt_ref and alt == axt_anc:
					flip = False
				# if the current allele is ancestral, the flipped count is for derived allele
				elif ref == axt_anc and alt == axt_ref:
					flip = True
				elif axt_ref == axt_anc and ref == axt_ref:
					flip = True
				elif axt_ref == axt_anc and alt == axt_ref:
					flip = False
				# what other cases could there be...?
				else:
					print(f'On chr{ch} position {pos}, ref, anc in axt: ({axt_ref}, {axt_anc}); in vcf: ({ref}, {alt}), set(axt_ref, axt_anc, ref, alt) = {set(list(axt_ref, axt_anc, ref, anc) )}.')
					sys.exit(1)
			# get genotype
			GT_id = _get_GT_index(l)
			x = 0; n = 0
			for i in individual_ids:
				GT = l[i].split(":")[GT_id]
				count = tuple(map( int, int_regex.findall(GT) ))
				# sanity check for bi-allelic
				try:
					assert set(count).issubset({0,1})
				except Exception as e:
					print("Warning: This script only applies to diploids.")
					print(e)
					sys.exit(1)
				x += sum(count)
				n += len(count)
			if flip:
				# remove x==n
				if x < n:
					# write output
					out.write(f"{pos}\t{pos * rec_rate}\t{n-x}\t{n}\n")
			else:
				# remove x==0
				if x > 0:
					# write output
					out.write(f"{pos}\t{pos * rec_rate}\t{x}\t{n}\n")
		
		# update axt_pos too
		while axt_pos <= pos:
			axt_pos_i += 1; axt_pos = aligned_pos[axt_pos_i]

	# end of file. Close both input and output
	vcf.close()
	out.close()



'''Read through the sequence alignment, polarize the frequency, and get genetic position from matching rec map'''
def parse_axt_and_vcf_w_recmap(ch, axtfile, vcffile, rec_map, rec_rate, outfile, pop_list):
	# load axt file
	print(f'Loading the alignment for chr{ch}...')
	axt_lib, aligned_pos = _load_alignment(axtfile, ch)

	# read and parse vcf:
	# write header of output
	out = open(outfile, 'w')
	out.write(f'position\tgenPos\tx\tn\n')

	# open vcf
	vcf = _open_large_file(vcffile, '.vcf')

	# open recomb map
	recmap = _open_large_file(rec_map)
	## start reading map:
	rec_l =  recmap.readline().strip().split("\t")
	# skip header if there's any
	if int_regex.findall(rec_l[1]) != rec_l[1]:
		print('Skipping header.')
		rec_l =  recmap.readline().strip().split("\t")
	rec_ch = rec_l[0]
	assert rec_ch in {ch, 'chr'+ch}
	r_pos, rate, cum_cM = map(float, rec_l[1:])
	if cum_cM == 0:
		rate = 0
	elif rate == 'NA':
		rate = rec_rate * 1e6 #typically in unit of cM/Mb
	last_r_pos, last_r_rate, last_cum_CM = 0, rate, 0
	last_pos, last_gen_pos = 0, 0

	# start reading vcf
	l = vcf.readline()
	## go to the header
	while l.startswith('##'):
		l = vcf.readline()
		continue
	## now l is vcf header
	## 0CHROM	1POS	2ID	3REF	4ALT	5QUAL	6FILTER	7INFO	8FORMAT	9:<individuals>
	individual_ids = _get_ids_from_vcf_header(l, pop_list)
	sampleSize = 2* len(individual_ids)
	try:
		#initiate indexing on aligned pos list
		axt_pos_i = 0 ; axt_pos = aligned_pos[axt_pos_i] ; axt_eof = False
	except Exception as e:
		print(e)
		print(aligned_pos, len(axt_lib))
		sys.exit(1)

	# start reading data
	while l != ['']:
		l = vcf.readline()
		if l == '':
			break
		l = l.strip().split('\t')
		# sanity check: it's the same chromosome
		try:
			assert l[0] in {ch, 'chr'+ch}
		except Exception as e:
			print(e)
			print(l[0], ch, {ch, 'chr'+ch})
			print(l[:11])
			sys.exit(1) 
		# check if pos in axt_pos & update axt_pos_i
		pos = int(l[1])
		while axt_pos < pos:
			axt_ref, axt_anc = axt_lib[axt_pos]
			# assume axt_ref; record as x=n=sampleSize
			if axt_ref != axt_anc:
				# find genetic position for axt_pos (in cM)
				# get genetic position (in cM)
				if axt_pos <= r_pos:
					# use the rec rate in last seg in recmap
					gen_pos = last_gen_pos + (pos - last_pos) * (last_r_rate)/1e6
				# otherwise, read another line
				else:
					while axt_pos > r_pos:
						# update and read another line from recmap
						last_r_pos, last_r_rate, last_cum_CM = r_pos, rate, cum_cM
						rec_l =  recmap.readline().strip().split("\t")
						if rec_l == ['']: #presumably end-of-file
							print(f'End-of-file for the recombination map. Assume uniform rate of {rec_rate} after {last_r_pos}.')
							r_pos, rate, cum_cM = 8e8, rec_rate*1e6 , (last_cum_CM + rec_rate*(8e8-last_r_pos))
							break
						rec_ch = rec_l[0]
						assert int_regex.findall(rec_ch)[0] == int_regex.findall(l[0])[0]
						r_pos, rate, cum_cM = map(float, rec_l[1:])
						if rate == 'NA':
							rate = rec_rate * 1e6 #typically in unit of cM/Mb
					if last_pos < last_r_pos:
						gen_pos = last_cum_CM + (axt_pos - last_r_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
					else:
						gen_pos = last_gen_pos + (axt_pos - last_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
				# update gen_pos
				last_pos, last_gen_pos = axt_pos, gen_pos
				# write output
				out.write(f"{axt_pos}\t{gen_pos}\t{sampleSize}\t{sampleSize}\n")
			# update axt_pos
			axt_pos_i += 1
			try:
				axt_pos = aligned_pos[axt_pos_i]
			except IndexError:
				print(f'No more recorded positions in axt.')
				axt_eof = True
				break
		# end of while-loop for axt
		if axt_eof:
			break
		## Now axt_pos >= pos
		# only process bi-allelic SNPs
		if l[3] in SNP and l[4] in SNP and "PASS" in l[6]:
			# skip if not mapped
			if pos not in axt_lib:
				continue
			# skip if not bi-allelic
			elif len(set([axt_lib[pos], l[3], l[4]] )) > 2:
				continue
			# polarize
			else:
				assert axt_pos == pos
				axt_ref, axt_anc = axt_lib[pos]
				ref, alt = l[3], l[4]
				# if the alternate allele is ancestral, current count is for derived allele
				if ref == axt_ref and alt == axt_anc:
					flip = False
				# if the current allele is ancestral, the flipped count is for derived allele
				elif ref == axt_anc and alt == axt_ref:
					flip = True
				elif axt_ref == axt_anc and ref == axt_ref:
					flip = True
				elif axt_ref == axt_anc and alt == axt_ref:
					flip = False
				# what other cases could there be...?
				else:
					print(f'On chr{ch} position {pos}, ref, anc in axt: ({axt_ref}, {axt_anc}); in vcf: ({ref}, {alt}), set(axt_ref, axt_anc, ref, alt) = {set(list(axt_ref, axt_anc, ref, anc) )}.')
					sys.exit(1)
			# get genotype
			## from "FORMAT" get how to obtain GT
			GT_id = _get_GT_index(l)
			x = 0; n = 0
			for i in individual_ids:
				GT = l[i].split(":")[GT_id]
				count = tuple(map( int, int_regex.findall(GT) ))
				# sanity check for bi-allelic
				try:
					assert set(count).issubset({0,1})
				except Exception as e:
					print("Warning: This script only applies to diploids.")
					print(e)
					sys.exit(1)
				x += sum(count)
				n += len(count)

			# get genetic position (in cM)
			if pos <= r_pos:
				# use the rec rate in last seg in recmap
				gen_pos = last_gen_pos + (pos - last_pos) * (last_r_rate)/1e6
			# otherwise, read another line of recmap
			else:
				while pos > r_pos:
					# update and read another line from recmap
					last_r_pos, last_r_rate, last_cum_CM = r_pos, rate, cum_cM
					rec_l =  recmap.readline().strip().split("\t")
					if rec_l == ['']: #presumably end-of-file
						print(f'End-of-file for the recombination map. Assume uniform rate of {rec_rate} after {last_r_pos}.')
						r_pos, rate, cum_cM = 8e8, rec_rate*1e6 , (last_cum_CM + rec_rate*(8e8-last_r_pos))
						break
					rec_ch = rec_l[0]
					assert int_regex.findall(rec_ch)[0] == int_regex.findall(l[0])[0]
					r_pos, rate, cum_cM = map(float, rec_l[1:])
					if rate == 'NA':
						rate = rec_rate * 1e6 #typically in unit of cM/Mb
				if last_pos < last_r_pos:
					gen_pos = last_cum_CM + (pos - last_r_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
				else:
					gen_pos = last_gen_pos + (pos - last_pos) / (r_pos - last_r_pos) * (cum_cM - last_cum_CM)
			# update gen_pos
			last_pos, last_gen_pos = pos, gen_pos

			if flip:
				# remove x==n
				if x < n:
					# write output
					out.write(f"{pos}\t{gen_pos}\t{n-x}\t{n}\n")
			else:
				# remove x==0
				if x > 0:
					# write output
					out.write(f"{pos}\t{gen_pos}\t{x}\t{n}\n")
		
		# update axt_pos too
		while axt_pos <= pos and axt_pos_i < len(aligned_pos):
			axt_pos_i += 1; axt_pos = aligned_pos[axt_pos_i]
		#read next line
		#l = vcf.readline()
		#print(l)
	# end of while loop through vcf
	# end of file. Close both input and output
	vcf.close()
	out.close()
	recmap.close()

#===================main body below===============
def main():
	import argparse
	parser = argparse.ArgumentParser()

	# file paths
	parser.add_argument('--vcf', dest = 'vcffile', help = 'Path and name of the vcf file. Format can be either .vcf or .vcf.gz.\n', required = True)
	parser.add_argument('-c', '--chr', dest = 'ch', help = 'ID of the chromosome. E.g. 2a for chr2a, 12 for chr12, etc.\n', required = True)
	parser.add_argument('-o','--output', dest = 'outfile', help = 'Path and name of the output file.', required = True)
	parser.add_argument('--ID_list', dest = 'pop_list', default = None, help = 'Path and name to the file containing list of sample IDs (identical to their column names in vcf) to be counted, separated by comma. If not provided, all samples in the vcf will be counted.\n')
	parser.add_argument('--axt', dest = 'axtfile', default = None, help = 'Path and name of the sequence alignment (in .axt or .axt.gz) for calling substitution and polarizing the allele frequency. If not provided, then the output will only be applicable to B_0maf.')
	parser.add_argument('--rec', dest = 'rec_rate', default = 1e-6, help = 'Recombination rate in cM/nt. Default value is 1e-6 cM/nt.')
	parser.add_argument('--rec_map', dest = 'rec_map', default = None, help = 'Path and name of the recombination map (hapmap format) of the same sequence. If not provided, a uniform recombination rate will be applied with a default rate of 1e-6 cM/nt. Use \"--rec\" to specify another rate.')

	
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	opt = parser.parse_args(sys.argv[1:])

	import time
	last = time.time()
	if opt.rec_map is not None:
		if opt.axtfile is not None:
			print(time.ctime(), f'Parsing BalLeRMix input for B_2 with genetic positions matching the recombination map. Assume genomic regions not covered by the map to recombine at {opt.rec_rate} cM/nt.')
			parse_axt_and_vcf_w_recmap(opt.ch, opt.axtfile, opt.vcffile, opt.rec_map, opt.rec_rate, opt.outfile, opt.pop_list)
		else:
			print(time.ctime(), f'Parsing BalLeRMix input for B_0maf with genetic positions matching the recombination map. Assume genomic region not covered by the map to recombine at {opt.rec_rate} cM/nt.')
			parse_vcf_w_recmap(opt.ch, opt.vcffile, opt.rec_map, opt.rec_rate, opt.outfile, opt.pop_list)
	else:
		if opt.axtfile is not None:
			print(time.ctime(), f'Parsing BalLeRMix input for B_2 with a uniform recombination rate of {opt.rec_rate} cM/nt...')
			parse_axt_and_vcf(opt.ch, opt.axtfile, opt.vcffile, opt.rec_rate, opt.outfile, opt.pop_list)
		else:
			print(time.ctime(), f'Parsing BalLeRMix input for B_0maf with a uniform recombination rate of {opt.rec_rate} cM/nt...')
			parse_vcf_only(opt.ch, opt.vcffile, opt.rec_rate, opt.outfile, opt.pop_list)

	print(time.ctime(),f'Parsing completed. Output file: {opt.outfile} Total time {time.time() - last}sec.')

if __name__ == '__main__':
	main()

