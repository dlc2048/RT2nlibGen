
import os
from src.endf.endf import ENDF

target_list = os.listdir('endf_viii0')
output_list = []
for file in target_list:
    print(file)
    file      = os.path.join('endf_viii0', file)
    endf_data = ENDF(file, verbose=False, read_only_header=False)
    ad = endf_data[2].ad
    if ad is not None:
        a = ad.getDistribution(19.95e6)
        if not a.isLegendre():
            print('tuple')

# initialize subprocess
print(1)
