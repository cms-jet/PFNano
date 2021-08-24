import FWCore.ParameterSet.Config as cms

pfMassDecorrelatedParticleNetPreprocessParamsAK15 = cms.PSet(
    input_names = cms.vstring(
        'pf_points',
        'pf_features',
        'pf_mask',
        'sv_points',
        'sv_features',
        'sv_mask'
    ),
    pf_features = cms.PSet(
        var_infos = cms.PSet(
            pfcand_VTX_ass = cms.PSet(
                median = cms.double(4),
                norm_factor = cms.double(0.3), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_abseta = cms.PSet(
                median = cms.double(0.7),
                norm_factor = cms.double(1.4), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_btagEtaRel = cms.PSet(
                median = cms.double(2),
                norm_factor = cms.double(0.4),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-5),
                upper_bound = cms.double(15),
            ),
            pfcand_btagJetDistVal = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(145),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-20),
                upper_bound = cms.double(1),
            ),
            pfcand_btagPParRatio = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-10),
                upper_bound = cms.double(100),
            ),
            pfcand_btagPtRatio = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(10),
            ),
            pfcand_btagSip3dSig = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(4e4),
            ),
            pfcand_btagSip3dVal = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(230),
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(1e5),
            ),
            pfcand_charge = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_dxy = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(250), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_dxysig = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1.6), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-2000),
                upper_bound = cms.double(2000),
            ),
            pfcand_dz = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(170), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_dzsig = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1.2), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-2000),
                upper_bound = cms.double(2000),
            ),
            pfcand_e_log_nopuppi = cms.PSet(
                median = cms.double(1),
                norm_factor = cms.double(0.6), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_etarel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_lostInnerHits = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_normchi2 = cms.PSet(
                median = cms.double(5),
                norm_factor = cms.double(0.2), 
                replace_inf_value = cms.double(300),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(300),
            ),
            pfcand_phirel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_pt_log_nopuppi = cms.PSet(
                median = cms.double(0.6),
                norm_factor = cms.double(0.6), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_quality = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(0.2), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(24),
        max_length = cms.uint32(150),
        var_names = cms.vstring(
            'pfcand_pt_log_nopuppi',
            'pfcand_e_log_nopuppi',
            'pfcand_etarel',
            'pfcand_phirel',
            'pfcand_abseta',
            'pfcand_charge',
            'pfcand_VTX_ass',
            'pfcand_lostInnerHits',
            'pfcand_normchi2',
            'pfcand_quality',
            'pfcand_dz',
            'pfcand_dzsig',
            'pfcand_dxy',
            'pfcand_dxysig',
            'pfcand_btagEtaRel',
            'pfcand_btagPtRatio',
            'pfcand_btagPParRatio',
            'pfcand_btagSip3dVal',
            'pfcand_btagSip3dSig',
            'pfcand_btagJetDistVal'
        ), 
    ),
    pf_mask = cms.PSet(
        var_infos = cms.PSet(
            pfcand_mask = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(24),
        max_length = cms.uint32(150),
        var_names = cms.vstring('pfcand_mask'), 
    ),
    pf_points = cms.PSet(
        var_infos = cms.PSet(
            pfcand_etarel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            pfcand_phirel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(24),
        max_length = cms.uint32(150),
        var_names = cms.vstring(
            'pfcand_etarel',
            'pfcand_phirel'
        )
    ),
    sv_features = cms.PSet(
        var_infos = cms.PSet(
            sv_abseta = cms.PSet(
                median = cms.double(0.7),
                norm_factor = cms.double(1.4), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_d3d = cms.PSet(
                median = cms.double(0.4),
                norm_factor = cms.double(0.3), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_d3dsig = cms.PSet(
                median = cms.double(5),
                norm_factor = cms.double(0.2), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(800),
            ),
            sv_dxy = cms.PSet(
                median = cms.double(0.25),
                norm_factor = cms.double(0.4), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_dxysig = cms.PSet(
                median = cms.double(5),
                norm_factor = cms.double(0.2), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1),
                upper_bound = cms.double(800),
            ),
            sv_etarel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_mass = cms.PSet(
                median = cms.double(3),
                norm_factor = cms.double(0.3), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_normchi2 = cms.PSet(
                median = cms.double(1),
                norm_factor = cms.double(0.5), 
                replace_inf_value = cms.double(1000),
                lower_bound = cms.double(-1000),
                upper_bound = cms.double(1000),
            ),
            sv_ntracks = cms.PSet(
                median = cms.double(3.0),
                norm_factor = cms.double(1.0), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_phirel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_pt_log = cms.PSet(
                median = cms.double(3),
                norm_factor = cms.double(0.6), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(1),
        max_length = cms.uint32(10),
        var_names = cms.vstring(
            'sv_pt_log',
            'sv_mass',
            'sv_etarel',
            'sv_phirel',
            'sv_abseta',
            'sv_ntracks',
            'sv_normchi2',
            'sv_dxy',
            'sv_dxysig',
            'sv_d3d',
            'sv_d3dsig'
        )
    ),
    sv_mask = cms.PSet(
        var_infos = cms.PSet(
            sv_mask = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(1),
        max_length = cms.uint32(10),
        var_names = cms.vstring('sv_mask')
    ),
    sv_points = cms.PSet(
        var_infos = cms.PSet(
            sv_etarel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            ),
            sv_phirel = cms.PSet(
                median = cms.double(0),
                norm_factor = cms.double(1), 
                replace_inf_value = cms.double(0.0),
                lower_bound = cms.double(-1e32),
                upper_bound = cms.double(1e32),
            )
        ),
        min_length = cms.uint32(1),
        max_length = cms.uint32(10),
        var_names = cms.vstring(
            'sv_etarel',
            'sv_phirel'
        )
    )
)
