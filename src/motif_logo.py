from weblogolib import *

# See https://github.com/WebLogo/weblogo/blob/master/weblogolib/__init__.py

# To create a PNG logo in python code do something like this:
fin = open('cap.fa')
seqs = read_seq_data(fin)
data = LogoData.from_seqs(seqs)
options = LogoOptions()
options.title = "A Logo Title"
format = LogoFormat(data, options)
png = png_formatter(data, format)

