import FWCore.ParameterSet.Config as cms
import os
import fnmatch

process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
#process.MessageLogger.destinations = cms.untracked.vstring("debug_rootupleChib.txt")
#process.MessageLogger.cout = process.MessageLogger.cerr
#process.MessageLogger.cout.threshold = cms.untracked.string("INFO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# input dir
#dir = '/data1/gdujany/chiB/CMSSW_5_3_3/src/SimChiB/crab/Multicrab_ChiB_1P_1_2_UpsilonPt_10M_each/ChiB_1P_1_UpsilonPt_10M/res/'


#inputfiles = ['file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_1_1_5YM.root']


inputfiles = [
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_14_1_gt0.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_18_1_0CG.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_19_1_ZXu.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_1_1_5YM.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_22_1_lED.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_24_1_Xtl.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_25_1_SVQ.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_2_1_2Nb.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_33_1_eqB.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_34_1_Swa.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_35_1_dED.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_36_1_wcy.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_38_1_Hta.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_39_1_SwG.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_3_1_A4L.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_41_1_FIQ.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_43_1_EEH.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_44_1_618.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_49_1_QJH.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_4_1_I5m.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_55_1_BN0.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_57_1_9ZW.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_58_1_HJT.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_5_1_gvC.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_61_1_hJv.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_63_1_dsK.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_66_1_aRc.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_68_1_n74.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_6_1_QZm.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_76_1_08B.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_79_1_I3V.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_7_1_uvN.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_80_1_zD1.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_85_1_IJy.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_87_1_opi.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_89_1_iFm.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_8_1_xzA.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_91_1_V8n.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_92_1_tQl.root',
    'file:/data2/data/degano/ponia_development/CMSSW_5_3_3/src/test/chic2012B_withTrigger/res/sel_chic_93_1_r9w.root',
        
    ]

#inputfiles = []
#for file in os.listdir(dir):
#    if fnmatch.fnmatch(file, 'testdimuon*.root'):
#        filename = 'file:'+dir+str(file)
#        inputfiles.append(filename)

process.source = cms.Source("PoolSource",
    	fileNames = cms.untracked.vstring( inputfiles )
)
process.TFileService = cms.Service("TFileService", 
	fileName = cms.string("rootuplaTest.root"),
	closeFileFast = cms.untracked.bool(True)
)


process.load('Ponia.RootupleChib.rootuplechic_cfi')

process.p = cms.Path(process.rootuple)
