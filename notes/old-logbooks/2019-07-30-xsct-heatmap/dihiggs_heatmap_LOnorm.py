from hfile import get_array
import numpy as np
import argparse
from matplotlib import pyplot as plt
import seaborn as sns

def load_heatmap(filename):
    ybins=['0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0']
    
    dist2d=np.zeros((11,10))
    i=0
    for y in ybins:
        ydist = get_array(filename, 'yjjcut-%s-mjj'%y)
        for j in range(len(ydist)):
            # print([x for x in ydist[:j+1,2]],':')
            # print(ydist[j,:2],sum([x for x in ydist[j:,2]]))
            ydist[j,2] = sum([x for x in ydist[j:,2]])
        dist2d[i,:] = ydist[:,2]
        i+=1
    return dist2d

# parser = argparse.ArgumentParser(description='Create a heatmap from histograms.')
# parser.add_argument('file', action='store', default=None,
#                     help='A top file with the histograms.')
# parser.add_argument('--output', '-o', type=str, default=None,
#                     help='The output pdf file.')
# args = parser.parse_args()

fnfact_nnlo='dihiggs-heatmap-fact-nnlo.top'
fnfact_nlo='dihiggs-heatmap-fact-nlo.top'
fn_lo='dihiggs-heatmap-lo.top'
fnnonfact='dihiggs-heatmap-nonfact.top'
#fnnonfact='pwg-LO.top'
#yval   = [0.0, 5.5]
yval   = [-0.25, 5.25]
ytic   = [0, 1, 2, 3, 4, 5]
#xval   = [0.0, 1000.0]
xval   = [-50.0, 950.0]
xtic   = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

dist2d_nnlo = load_heatmap(fnfact_nnlo)
dist2d_nlo = load_heatmap(fnfact_nlo)
dist2d_lo = load_heatmap(fn_lo)
dist2d = dist2d_nnlo-dist2d_nlo
dist2d_nonfact = load_heatmap(fnnonfact)


fig = plt.figure(figsize=(14,6))
plt.rcParams.update({'font.size': 16})
plt.xlabel('$m_{jj,\mathrm{cut}}$ [GeV]')
plt.ylabel('$\Delta y_{jj,\mathrm{cut}}$')
plt.title('proVBFHH v1.2.0')
# im = plt.imshow(avg_jetimg.transpose(), origin='lower', aspect='auto',
#                 extent=xval+yval, cmap=plt.get_cmap('Reds'), vmin=0.0, vmax=0.02)
ratio=100*dist2d_nonfact/dist2d_lo
im = plt.imshow(ratio, origin='lower', aspect='auto',
                extent=xval+yval,cmap=plt.get_cmap('seismic'), vmin=-4.2, vmax=4.2)
cbar = fig.colorbar(im, label='[%]')

ysize = dist2d.shape[0]
xsize = dist2d.shape[1]
jump_x = (xval[1] - xval[0]) / (2.0 * xsize)
jump_y = (yval[1] - yval[0]) / (2.0 * ysize)
x_positions = np.linspace(start=xval[0], stop=xval[1], num=xsize, endpoint=False)
y_positions = np.linspace(start=yval[0], stop=yval[1], num=ysize, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = '%.2f'%ratio[y_index, x_index]
        text_x = x + jump_x
        text_y = y + jump_y
        plt.text(text_x, text_y, label, color='black', ha='center', va='center', fontsize=13)
plt.xticks(xtic)
plt.yticks(ytic)

plt.savefig('dihiggs_heatmap_nonfact_normLO_new.pdf', bbox_inches='tight')

fig = plt.figure(figsize=(14,6))
plt.rcParams.update({'font.size': 16})
plt.xlabel('$m_{jj,\mathrm{cut}}$ [GeV]')
plt.ylabel('$\Delta y_{jj,\mathrm{cut}}$')
plt.title('proVBFHH v1.2.0')
# im = plt.imshow(avg_jetimg.transpose(), origin='lower', aspect='auto',
#                 extent=xval+yval, cmap=plt.get_cmap('Reds'), vmin=0.0, vmax=0.02)
ratio=100*dist2d/dist2d_lo
im = plt.imshow(ratio, origin='lower', aspect='auto',
                extent=xval+yval,cmap=plt.get_cmap('seismic'), vmin=-4.2, vmax=4.2)
cbar = fig.colorbar(im, label='[%]')

ysize = dist2d.shape[0]
xsize = dist2d.shape[1]
jump_x = (xval[1] - xval[0]) / (2.0 * xsize)
jump_y = (yval[1] - yval[0]) / (2.0 * ysize)
x_positions = np.linspace(start=xval[0], stop=xval[1], num=xsize, endpoint=False)
y_positions = np.linspace(start=yval[0], stop=yval[1], num=ysize, endpoint=False)

for y_index, y in enumerate(y_positions):
    for x_index, x in enumerate(x_positions):
        label = '%.2f'%ratio[y_index, x_index]
        text_x = x + jump_x
        text_y = y + jump_y
        plt.text(text_x, text_y, label, color='white', ha='center', va='center', fontsize=13)
plt.xticks(xtic)
plt.yticks(ytic)

plt.savefig('dihiggs_heatmap_fact_normLO.pdf', bbox_inches='tight')
