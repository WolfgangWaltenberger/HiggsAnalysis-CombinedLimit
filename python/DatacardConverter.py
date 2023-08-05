import ROOT as r
import json, math
from HiggsAnalysis.CombinedLimit.DatacardParser import *

def checkBin(v, neg):
    
    if v < 1E-10 and neg: return 0.
    else: return v

def convertCard(cardName, f, opts, outName, bbl, normshape, prune, neg):

    card = {'channels': [], 'observations': [], 'measurements': [], 'version': '1.0.0'}
    
    zer = 1E-20

    with open(cardName, 'r') as f:
        dc = parseCard(f, opts)
        sig = dc.signals[0]
        meas = []
        for ich, chname in enumerate(dc.bins):
            ch = {'name': chname, 'samples': []}
            fdata = list(dc.shapeMap[chname].values())[0]
            fname = fdata[0]
            hnom = fdata[1]
            hsys = fdata[2]
            fr = r.TFile(fname, 'READ')
            for k in dc.keyline:
                if k[0] != chname: continue
                s = k[1]
                hnomName = hnom.replace('$PROCESS', s)
                h = fr.Get(hnomName)                
                hInteg = h.Integral()
                if hInteg < 1E-10: continue
                hData = [checkBin(h.GetBinContent(ib+1), neg) for ib in range(h.GetXaxis().GetNbins())]
                hNorm = h.Clone('hNorm')
                hNorm.Scale(1./hInteg)
                hNormData = [checkBin(hNorm.GetBinContent(ib+1), neg) for ib in range(hNorm.GetXaxis().GetNbins())]
                if h == None: continue
                ch['samples'].append({})
                data = []
                nBins = h.GetXaxis().GetNbins()
                for ib in range(nBins):
                    data.append(checkBin(h.GetBinContent(ib+1)+zer, neg))
                ch['samples'][-1]['name'] = s
                ch['samples'][-1]['data'] = data
                ch['samples'][-1]['modifiers'] = []
            
                if s == sig:
                    ch['samples'][-1]['modifiers'].append({'data': None, 'name': 'r_'+sig, 'type': 'normfactor'})

                for syst in dc.systs:
                    systname = syst[0]
                    systtype = syst[2]
                    systdict = syst[4]
                    if s not in systdict[chname].keys(): continue
                    systfact = systdict[chname][s]

                    if systfact != 0.0:

                        systIncl, systInclNorm = False, False
                        if systtype in ['lnN', 'lnU']:
                            if type(systfact) != list and abs(systfact-1.0) > 1E-10:
                                systdata = {'name': systname}
                                systdata['type'] = 'normsys'
                                normsysup = systfact
                                normsysdo = 1.0/systfact
#                                normsysdo = 1.0-abs(systfact-1.0)
                                systdata['data'] = {'hi': normsysup, 'lo': normsysdo}
                                systIncl = True
                            elif type(systfact) is list:
                                systdata = {'name': systname}
                                systdata['type'] = 'normsys'
                                normsysup = systfact[1]
                                normsysdo = systfact[0]
                                systdata['data'] = {'hi': normsysup, 'lo': normsysdo}
                                systIncl = True
                        elif systtype in ['shape']:
                            systdata = {'name': systname}
                            systdata['type'] = 'histosys'
                            hsysNameUp = hsys.replace('$PROCESS', s).replace('$SYSTEMATIC', systname+'Up')
                            hsysNameDown = hsys.replace('$PROCESS', s).replace('$SYSTEMATIC', systname+'Down')
                            hsysUp = fr.Get(hsysNameUp)
                            hsysDown = fr.Get(hsysNameDown)
                            hsysNormUp = hsysUp.Clone('hsysNormUp')
                            hsysNormDown = hsysDown.Clone('hsysNormDown')
                            hsysIntegUp = hsysUp.Integral()
                            hsysIntegDown = hsysDown.Integral()
                            hasNorm = bool(abs(hsysIntegUp-hsysIntegDown) > 1E-10)
                            normUp = hsysIntegUp/hInteg
                            normDown = hsysIntegDown/hInteg
                            if normshape or '_splitns' in systname:
                                hsysUp.Scale(1./normUp)
                                hsysDown.Scale(1./normDown)
                            hsysNormUp.Scale(1./hsysIntegUp)
                            hsysNormDown.Scale(1./hsysIntegDown)
                            hsysDataUp = [checkBin(hsysUp.GetBinContent(ib+1), neg) for ib in range(hsysUp.GetXaxis().GetNbins())]
                            hsysDataDown = [checkBin(hsysDown.GetBinContent(ib+1), neg) for ib in range(hsysDown.GetXaxis().GetNbins())]
                            hsysNormDataUp = [checkBin(hsysNormUp.GetBinContent(ib+1), neg) for ib in range(hsysNormUp.GetXaxis().GetNbins())]
                            hsysNormDataDown = [checkBin(hsysNormDown.GetBinContent(ib+1), neg) for ib in range(hsysNormDown.GetXaxis().GetNbins())]
                            diffShapeUp, diffShapeDown = 0., 0.
                            for ib in range(nBins):
                                nup = abs(hsysNormDataUp[ib])+abs(hNormData[ib])
                                vup = 2.*abs(hsysNormDataUp[ib]-hNormData[ib])
                                diffShapeUp += vup/nup if nup > 0 else 0.
                                ndown = abs(hsysNormDataDown[ib])+abs(hNormData[ib])
                                vdown = 2.*abs(hsysNormDataDown[ib]-hNormData[ib])
                                diffShapeDown += vdown/ndown if ndown > 0 else 0.
                            hasShape = bool(diffShapeUp > 1E-4 and diffShapeDown > 1E-4) or prune
                            isFake = (len(hsysDataUp) == 1 and abs(hsysDataUp[0]-hsysDataDown[0]) < 1E-4)
                            if systfact != 1.0:
                                print('Warning: an additional shape normalization factor found')
                                for ib in range(nBins):
                                    hsysDataUp[ib] = checkBin((hsysDataUp[ib]-hData[ib])*systfact+hData[ib], neg)
                                    hsysDataDown[ib] = checkBin((hsysDataDown[ib]-hData[ib])*systfact+hData[ib], neg)
                            if (hasNorm and normshape) or ('_splitns' in systname):
                                systdatanorm = {'name': systname}
                                if '_splitns' not in systname: systdatanorm['name'] += '_mergedns'
                                systdatanorm['type'] = 'normsys'
                                systdatanorm['data'] = {'hi': normUp, 'lo': normDown}
                                systInclNorm = True
#                                if abs(sum(hsysNormDataUp)-sum(hsysNormDataDown)) > 1E-7:
                                if True:
                                    if '_splitns' not in systname: systdata['name'] += '_mergedns'
                                    else:
                                        systdatanorm['name'] = systdatanorm['name'].replace('_splitns', '')
                                        if hasShape and not isFake:
                                            systdata['name'] = systdata['name'].replace('_splitns', '')
                                            systdata['data'] = {'hi_data': hsysDataUp, 'lo_data': hsysDataDown}
                                    systIncl = True
                            if hasShape:
                                systdata['data'] = {'hi_data': hsysDataUp, 'lo_data': hsysDataDown}
                                systIncl = True

                        if systIncl: 
                            if 'data' in systdata.keys():
                                ch['samples'][-1]['modifiers'].append(systdata)
                        if systInclNorm and abs(systdatanorm['data']['hi']-systdatanorm['data']['lo']) > 1E-4: ch['samples'][-1]['modifiers'].append(systdatanorm)

                if chname in dc.binParFlags.keys() and dc.binParFlags[chname][1] and dc.binParFlags[chname][0] >= 0:
                    if bbl:
                        systname = 'prop_bin'+chname
                        systdata = {'name': systname}
                        systdata['type'] = 'staterror'
                        hstat = [h.GetBinError(ib+1) for ib in range(h.GetXaxis().GetNbins())]
                        for iv in range(len(hstat)):
                            if data[iv] < 1E-10: hstat[iv] = 0.
                        systdata['data'] = hstat
                        ch['samples'][-1]['modifiers'].append(systdata)
                    else:
                        systname = 'prop_bin'+chname+'_'+s
                        systdata = {'name': systname}
                        systdata['type'] = 'shapesys'
                        hstat = [h.GetBinError(ib+1) for ib in range(h.GetXaxis().GetNbins())]
                        for iv in range(len(hstat)):
                            if data[iv] < 1E-10: hstat[iv] = 0.
                        systdata['data'] = hstat
                        ch['samples'][-1]['modifiers'].append(systdata)                        
                
                for rp in dc.rateParams.keys():
                    specs = rp.split('AND')
                    rpchan = specs[0]
                    rpproc = specs[1]
                    if chname == rpchan and s == rpproc:
                        rpname = dc.rateParams[rp][0][0][0]
                        ch['samples'][-1]['modifiers'].append({'data': None, 'name': rpname, 'type': 'normfactor'})
                        found = False
                        for m in meas:
                            if m['name'] == rpname:
                                found = True
                                break
                        if not found: meas.append({'name': rpname})
            
            hdata = fr.Get(hnom.replace('$PROCESS', 'data_obs'))
            hobs = [hdata.GetBinContent(ib+1) for ib in range(hdata.GetXaxis().GetNbins())]
            card['observations'].append({'name': chname, 'data': hobs})

            card['channels'].append(ch)

        par = {'bounds': [[-10.0, 10.0]], 'fixed': False, 'name': 'r_'+sig}
        card['measurements'] = [{'config': {'parameters': [par], 'poi': 'r_'+sig}, 'name': 'meas'}]
        for m in meas:
            par = {'bounds': [[-10.0, 10.0]], 'fixed': False, 'name': m['name'], 'inits': [1.0]}
            card['measurements'][0]['config']['parameters'].append(par)

    json.dump(card, open(outName+'.json', 'w'), indent=4)
