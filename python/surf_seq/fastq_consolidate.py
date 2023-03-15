import glob

# Define folder to use
#temp_dir = str(input('Directory of files'))

# Read in files
for i in range(96):
	i += 1
	if i < 10:
		barcode = 'barcode0' + str(i)
	else:
		barcode = 'barcode' + str(i)
	temp_dir = 'barcoding/' + barcode
	filenames = sorted(glob.glob(temp_dir+ '/*.fastq'))
	#print(filenames)

	# Open new file
	consol_file = open(temp_dir + '/' + barcode + '.fastq', 'w')
			
	# Iterate over files
	for temp_file in filenames:
		print(temp_file)
		temp_f = open(temp_file, 'r')
		temp_lines = temp_f.readlines()
		temp_f.close()
		for line in temp_lines:
			consol_file.write(line)
	consol_file.close()