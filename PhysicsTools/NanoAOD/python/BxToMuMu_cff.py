from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms
import re

# can use cms.PSet.clone() method instead
def merge_psets(*argv):
    result = cms.PSet()
    for pset in argv:
        if isinstance(pset, cms._Parameterizable):
            for name in pset.parameters_().keys():
                value = getattr(pset,name)
                type = value.pythonTypeName()
                setattr(result,name,value)
    return result

def copy_pset(pset,replace_dict):
    result = cms.PSet()
    for name in pset.parameters_().keys():
        new_name = name
        for pattern,repl in replace_dict.items():
            new_name = re.sub(pattern,repl,new_name)
        if getattr(pset, name).pythonTypeName() in ('cms.PSet', 'cms.untracked.PSet'):
            value = getattr(pset, name).clone()
            for sname in value.parameters_().keys():
                svalue = getattr(value,sname)
                stype = svalue.pythonTypeName()
                if stype in ('cms.string', 'cms.untracked.string'):
                    new_svalue = svalue.value()
                    for pattern,repl in replace_dict.items():
                        new_svalue = re.sub(pattern,repl,new_svalue)
                    setattr(value,sname,new_svalue)
            setattr(result,new_name,value)
    return result

BxToMuMu = cms.EDProducer("BxToMuMuProducer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection=cms.InputTag("packedPFCandidates"),
    prunedGenParticleCollection = cms.InputTag("prunedGenParticles"),
    MuonMinPt=cms.double(1.),
    MuonMaxEta=cms.double(2.4),
    KaonMinPt=cms.double(1.),
    KaonMaxEta=cms.double(2.4),
    KaonMinDCASig=cms.double(-1.),
    DiMuonChargeCheck=cms.bool(False),
    minBKmmMass = cms.double(4.5),
    maxBKmmMass = cms.double(6.0),
    maxTwoTrackDOCA = cms.double(0.1),
    isMC = cms.bool(False)
)

BxToMuMuMc = BxToMuMu.clone( isMC = cms.bool(True) ) 

kinematic_pset = cms.PSet(
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic vertex fit validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit vertex probability"),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic vertex refitted mass"),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic vertex refitted mass error"),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit vertex displacement in XY plane"),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit vertex displacement significance in XY plane"),
    kin_cosAlpha = Var("userFloat('kin_cosAlpha')",    float, doc = "Kinematic fit: cosAlpha"),
)

BxToMuMuDiMuonTableVariables = merge_psets(
    cms.PSet(
        mu1_index    = Var("userInt('mu1_index')",         int,   doc = "Index of corresponding leading muon"),
        mu2_index    = Var("userInt('mu2_index')",         int,   doc = "Index of corresponding subleading muon"),
        mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
        # Kalman Fit
        kal_valid    = Var("userInt('kalman_valid')",      int,   doc = "Kalman vertex fit validity"),
        kal_vtx_prob = Var("userFloat('kalman_vtx_prob')", float, doc = "Kalman fit vertex probability"),
        kal_mass     = Var("userFloat('kalman_mass')",     float, doc = "Kalman vertex refitted mass"),
        kal_lxy      = Var("userFloat('kalman_lxy')",      float, doc = "Kalman fit vertex displacement in XY plane"),
        kal_slxy     = Var("userFloat('kalman_sigLxy')",   float, doc = "Kalman fit vertex displacement significance in XY plane"),
        ),
    kinematic_pset
)

BxToMuMuDiMuonMcTableVariables = merge_psets(
    BxToMuMuDiMuonTableVariables,
    cms.PSet(
        gen_mu1_pdgId  = Var("userInt('gen_mu1_pdgId')",    int,   doc = "Gen match: first muon pdg Id"),
        gen_mu1_mpdgId = Var("userInt('gen_mu1_mpdgId')",   int,   doc = "Gen match: first muon mother pdg Id"),
        gen_mu1_pt     = Var("userFloat('gen_mu1_pt')",     float, doc = "Gen match: first muon pt"),
        gen_mu2_pdgId  = Var("userInt('gen_mu2_pdgId')",    int,   doc = "Gen match: second muon pdg Id"),
        gen_mu2_mpdgId = Var("userInt('gen_mu2_mpdgId')",   int,   doc = "Gen match: second muon mother pdg Id"),
        gen_mu2_pt     = Var("userFloat('gen_mu2_pt')",     float, doc = "Gen match: second muon pt"),
        gen_pdgId      = Var("userInt('gen_pdgId')",        int,   doc = "Gen match: dimuon pdg Id"),
        gen_mpdgId     = Var("userInt('gen_mpdgId')",       int,   doc = "Gen match: dimuon mother pdg Id"),
        gen_mass       = Var("userFloat('gen_mass')",       float, doc = "Gen match: dimuon mass"),
        gen_pt         = Var("userFloat('gen_pt')",         float, doc = "Gen match: dimuon pt"),
        ),
)

BxToMuMuDiMuonTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMu","DiMuon"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuDiMuonTableVariables
)

BxToMuMuDiMuonMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuMc","DiMuon"),
    cut=cms.string(""),
    name=cms.string("mm"),
    doc=cms.string("Dimuon Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuDiMuonMcTableVariables
)

BxToMuMuBToKmumuTableVariables =  merge_psets(
    copy_pset(kinematic_pset,{"kin_":"nomc_"}),
    copy_pset(kinematic_pset,{"kin_":"jpsimc_"}),
    cms.PSet(
        mm_index = Var("userInt('mm_index')",                 int,   doc = "Index of dimuon pair"),
        kaon_charge = Var("userInt('kaon_charge')",           int,   doc = "kaon charge"),
        kaon_pt     = Var("userFloat('kaon_pt')",             float, doc = "kaon pt"),
        kaon_eta    = Var("userFloat('kaon_eta')",            float, doc = "kaon eta"),
        kaon_phi    = Var("userFloat('kaon_phi')",            float, doc = "kaon phi")
    )
)

BxToMuMuBToKmumuMcTableVariables = merge_psets(
    BxToMuMuBToKmumuTableVariables,
    cms.PSet(
        gen_kaon_pdgId  = Var("userInt('gen_kaon_pdgId')",    int,   doc = "Gen match: kaon pdg Id"),
        gen_kaon_mpdgId = Var("userInt('gen_kaon_mpdgId')",   int,   doc = "Gen match: kaon mother pdg Id"),
        gen_kaon_pt     = Var("userFloat('gen_kaon_pt')",     float, doc = "Gen match: kaon pt"),
        gen_pdgId       = Var("userInt('gen_pdgId')",         int,   doc = "Gen match: kmm pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float, doc = "Gen match: kmm mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float, doc = "Gen match: kmm pt"),
    )
)
        

BxToMuMuBToKmumuTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMu","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKmumuTableVariables
)

BxToMuMuBToKmumuMcTable = cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("BxToMuMuMc","BToKmumu"),
    cut=cms.string(""),
    name=cms.string("bkmm"),
    doc=cms.string("BToKmumu Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = BxToMuMuBToKmumuMcTableVariables
)

BxToMuMuSequence   = cms.Sequence(BxToMuMu)
BxToMuMuMcSequence = cms.Sequence(BxToMuMuMc)
BxToMuMuTables     = cms.Sequence(BxToMuMuDiMuonTable   * BxToMuMuBToKmumuTable)
BxToMuMuMcTables   = cms.Sequence(BxToMuMuDiMuonMcTable * BxToMuMuBToKmumuMcTable)
