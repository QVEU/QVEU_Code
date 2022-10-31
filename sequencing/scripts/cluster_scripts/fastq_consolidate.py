import glob

# Define folder to use
#temp_dir = str(input('Directory of files'))

# Read in files
barcode = str(input('barcode##'))
temp_dir = '../data/2022-01-22_nanopore_data/barcode' + barcode
filenames = sorted(glob.glob(temp_dir+ '/*.fastq'))
print()
print(filenames)

# Open new file
consol_file = open(temp_dir + '/barcode' + barcode + '.fastq', 'w')

# Iterate over files
for temp_file in filenames:
	print(temp_file)
	temp_f = open(temp_file, 'r')
	temp_lines = temp_f.readlines()
	temp_f.close()
	for line in temp_lines:
		consol_file.write(line)
consol_file.close()