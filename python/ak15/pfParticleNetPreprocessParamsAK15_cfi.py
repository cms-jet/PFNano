import FWCore.ParameterSet.Config as cms

pfParticleNetPreprocessParamsAK15 = cms.PSet(
    input_names = cms.vstring(
        'pf_points', 
        'pf_features', 
        'pf_mask', 
        'sv_points', 
        'sv_features', 
        'sv_mask'
    ),
    pf_features = cms.PSet(
        input_shape = cms.vuint32(1, 25, 100),
        var_infos = cms.PSet(
            pfcand_VTX_ass = cms.PSet(
                median = cms.double(7.0),
                norm_factor = cms.double(0.5)
            ),
            pfcand_abseta = cms.PSet(
                median = cms.double(0.704061985016),
                norm_factor = cms.double(1.42180856942)
            ),
            pfcand_btagEtaRel = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(0.426671793464)
            ),
            pfcand_btagJetDistVal = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(145.473299395)
            ),
            pfcand_btagPParRatio = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.01681132268)
            ),
            pfcand_btagPtRatio = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.9316753054)
            ),
            pfcand_btagSip3dSig = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(0.964987250714)
            ),
            pfcand_btagSip3dVal = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(230.394748655)
            ),
            pfcand_charge = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_dxy = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(258.422706714)
            ),
            pfcand_dxysig = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.86056529822)
            ),
            pfcand_dz = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(171.380755873)
            ),
            pfcand_dzsig = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.38224973484)
            ),
            pfcand_e_log_nopuppi = cms.PSet(
                median = cms.double(0.960983216763),
                norm_factor = cms.double(0.602489542217)
            ),
            pfcand_etarel = cms.PSet(
                median = cms.double(-0.0385179147124),
                norm_factor = cms.double(2.00600043717)
            ),
            pfcand_isChargedHad = cms.PSet(
                median = cms.double(1.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_isEl = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_isGamma = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_isMu = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_isNeutralHad = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_lostInnerHits = cms.PSet(
                median = cms.double(-1.0),
                norm_factor = cms.double(1.0)
            ),
            pfcand_normchi2 = cms.PSet(
                median = cms.double(999.0),
                norm_factor = cms.double(0.001001001001)
            ),
            pfcand_phirel = cms.PSet(
                median = cms.double(-0.00015363653074),
                norm_factor = cms.double(2.07833165902)
            ),
            pfcand_pt_log_nopuppi = cms.PSet(
                median = cms.double(0.577863454819),
                norm_factor = cms.double(0.614248987126)
            ),
            pfcand_quality = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(0.2)
            )
        ),
        var_length = cms.uint32(100),
        var_names = cms.vstring(
            'pfcand_pt_log_nopuppi', 
            'pfcand_e_log_nopuppi', 
            'pfcand_etarel', 
            'pfcand_phirel', 
            'pfcand_abseta', 
            'pfcand_charge', 
            'pfcand_isMu', 
            'pfcand_isEl', 
            'pfcand_isChargedHad', 
            'pfcand_isGamma', 
            'pfcand_isNeutralHad', 
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
        replace_inf_value = cms.double(0.0),
    ),
    pf_mask = cms.PSet(
        input_shape = cms.vuint32(1, 1, 100),
        var_infos = cms.PSet(
            pfcand_mask = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            )
        ),
        var_length = cms.uint32(100),
        var_names = cms.vstring('pfcand_mask')
    ),
    pf_points = cms.PSet(
        input_shape = cms.vuint32(1, 2, 100),
        var_infos = cms.PSet(
            pfcand_etarel = cms.PSet(
                median = cms.double(-0.0385179147124),
                norm_factor = cms.double(2.00600043717)
            ),
            pfcand_phirel = cms.PSet(
                median = cms.double(-0.00015363653074),
                norm_factor = cms.double(2.07833165902)
            )
        ),
        var_length = cms.uint32(100),
        var_names = cms.vstring(
            'pfcand_etarel', 
            'pfcand_phirel'
        )
    ),
    sv_features = cms.PSet(
        input_shape = cms.vuint32(1, 12, 7),
        var_infos = cms.PSet(
            sv_abseta = cms.PSet(
                median = cms.double(0.692841857672),
                norm_factor = cms.double(1.36117543642)
            ),
            sv_costhetasvpv = cms.PSet(
                median = cms.double(0.999670803547),
                norm_factor = cms.double(170.088769034)
            ),
            sv_d3d = cms.PSet(
                median = cms.double(0.378370791674),
                norm_factor = cms.double(0.351689099492)
            ),
            sv_d3dsig = cms.PSet(
                median = cms.double(6.71155357361),
                norm_factor = cms.double(0.0288290198784)
            ),
            sv_dxy = cms.PSet(
                median = cms.double(0.257832497358),
                norm_factor = cms.double(0.445355413049)
            ),
            sv_dxysig = cms.PSet(
                median = cms.double(6.6999809742),
                norm_factor = cms.double(0.0288460696637)
            ),
            sv_etarel = cms.PSet(
                median = cms.double(-0.019983493723),
                norm_factor = cms.double(3.68221968779)
            ),
            sv_mass = cms.PSet(
                median = cms.double(1.18094140291),
                norm_factor = cms.double(0.531873999853)
            ),
            sv_normchi2 = cms.PSet(
                median = cms.double(0.756912440062),
                norm_factor = cms.double(0.741447745478)
            ),
            sv_ntracks = cms.PSet(
                median = cms.double(3.0),
                norm_factor = cms.double(1.0)
            ),
            sv_phirel = cms.PSet(
                median = cms.double(0.000604838802246),
                norm_factor = cms.double(3.76744890394)
            ),
            sv_pt_log = cms.PSet(
                median = cms.double(3.26131427288),
                norm_factor = cms.double(0.69510485069)
            )
        ),
        var_length = cms.uint32(7),
        var_names = cms.vstring(
            'sv_pt_log', 
            'sv_mass', 
            'sv_phirel', 
            'sv_etarel', 
            'sv_abseta', 
            'sv_ntracks', 
            'sv_normchi2', 
            'sv_dxy', 
            'sv_dxysig', 
            'sv_d3d', 
            'sv_d3dsig', 
            'sv_costhetasvpv'
        )
    ),
    sv_mask = cms.PSet(
        input_shape = cms.vuint32(1, 1, 7),
        var_infos = cms.PSet(
            sv_mask = cms.PSet(
                median = cms.double(0.0),
                norm_factor = cms.double(1.0)
            )
        ),
        var_length = cms.uint32(7),
        var_names = cms.vstring('sv_mask')
    ),
    sv_points = cms.PSet(
        input_shape = cms.vuint32(1, 2, 7),
        var_infos = cms.PSet(
            sv_etarel = cms.PSet(
                median = cms.double(-0.019983493723),
                norm_factor = cms.double(3.68221968779)
            ),
            sv_phirel = cms.PSet(
                median = cms.double(0.000604838802246),
                norm_factor = cms.double(3.76744890394)
            )
        ),
        var_length = cms.uint32(7),
        var_names = cms.vstring(
            'sv_phirel', 
            'sv_etarel'
        )
    )
)