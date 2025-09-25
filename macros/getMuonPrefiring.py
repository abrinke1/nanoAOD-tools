import numpy as np
import ROOT as R
in_file = R.TFile('data/prefire_maps/L1MuonPrefiringParametriations.root','R')
eras = ['2016preVFP','2016postVFP','20172018']
edg = ['0.0','0.2','0.3','0.55','0.83','1.24','1.4','1.6','1.8','2.1','2.25','2.4']
for era in eras:
    for i in range(1,len(edg)):
        func = in_file.Get('L1prefiring_muonparam_%sTo%s_%s' % (edg[i-1], edg[i], era))
        p0 = func.GetParameter(0)
        p1 = func.GetParameter(1)
        p2 = func.GetParameter(2)
        print('In %sTo%s_%s, wgt = 1 - (%.4f / (exp((x - %.1f) / %.2f) + 1))' % (edg[i-1], edg[i], era, p2, p0, p1))
        for pt in [10,25,1000]:
            print('  * pT = %d --> wgt = %.4f' % (pt, 1 - (p2 / (np.exp((pt - p0) / p1) + 1))))
