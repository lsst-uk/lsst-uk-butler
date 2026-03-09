import getopt
import os
import sys
import glob

from astropy.table import Table
import pandas as pd
from collections import defaultdict
from dotenv import load_dotenv

#-----------------------------------------------------------------------------

load_dotenv(os.path.join(os.environ['HOME'], '.env'))
defDataDir = 'dp1/LSSTComCam'
baseDir = os.path.join(os.environ["DATA"])
csvBaseDir = '../csv_files'
ecsvBaseDir = '../ecsv'
groups = {"deep_coadd": ['filename', 'band', 'skymap', 'tract', 'patch'],
          "visit_image": ['filename', 'band', 'instrument', 'day_obs',
                          'detector', 'physical_filter', 'visit'],
          "difference_image": ['filename', 'band', 'instrument', 'day_obs',
                               'detector', 'physical_filter', 'visit'],
          "template_coadd": ['filename', 'band', 'skymap', 'tract', 'patch'],
          "raw": ['filename', 'band', 'instrument', 'day_obs', 'detector',
                  'group', 'physical_filter', 'exposure'],
          "object": ['filename', 'skymap', 'tract'],
          "source": ['filename', 'band', 'instrument', 'day_obs',
                     'physical_filter', 'visit'],
          "object_forced_source": ['filename', 'skymap', 'tract', 'patch'],
          "dia_object": ['filename', 'skymap', 'tract'],
          "dia_source": ['filename', 'skymap', 'tract'],
          "dia_object_forced_source": ['filename', 'skymap', 'tract', 'patch'],
          "ss_object": ['filename',],
          "ss_source": ['filename',],
          "visit_table": ['filename', 'instrument'],
          "visit_detector_table": ['filename', 'instrument'],
          "deep_coadd_background": ['filename', 'band', 'skymap', 'tract', 'patch'],
          "deep_coadd_n_image": ['filename', 'band', 'skymap', 'tract', 'patch'],
          "visit_summary": ['filename', 'band', 'instrument', 'day_obs', 'physical_filter', 'visit'],
          "visit_image_background": ['filename', 'band', 'instrument', 'day_obs', 'detector', 'physical_filter', 'visit'],
          "the_monster_20250219": ['filename', 'htm7'],
          "deepCoadd_dcr_ddec_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_dcr_dra_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_dcr_e1_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_dcr_e2_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_epoch_consolidated_map_max": ['filename', 'band', 'skymap'],
          "deepCoadd_epoch_consolidated_map_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_epoch_consolidated_map_min": ['filename', 'band', 'skymap'],
          "deepCoadd_exposure_time_consolidated_map_sum": ['filename', 'band', 'skymap'],
          "deepCoadd_psf_e1_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_psf_e2_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_psf_maglim_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_psf_size_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_sky_background_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "deepCoadd_sky_noise_consolidated_map_weighted_mean": ['filename', 'band', 'skymap'],
          "mergeObjectMeasurement_metadata": ['filename', 'skymap', 'tract', 'patch'],
          "refCatObjectTract_metadata": ['filename', 'skymap', 'tract']
          }
excludedDT = {
              "analyzeDiaSourceTableTract_gatherResourceUsage_log": "json",
              "analyzeDiaSourceTableTract_gatherResourceUsage_metadata": "json",
              "analyzeDiaSourceTableTract_log": "json",
              "analyzeDiaSourceTableTract_metadata": "json",
              "analyzeObjectTableCore_gatherResourceUsage_log": "json",
              "analyzeObjectTableCore_gatherResourceUsage_metadata": "json",
              "analyzeObjectTableCore_log": "json",
              "analyzeObjectTableCore_metadata": "json",
              "analyzeObjectTableSurveyCore_gatherResourceUsage_log": "json",
              "analyzeObjectTableSurveyCore_gatherResourceUsage_metadata": "json",
              "analyzeObjectTableSurveyCore_log": "json",
              "analyzeObjectTableSurveyCore_metadata": "json",
              "analyzeRecalibratedStarAssociation_gatherResourceUsage_log": "json",
              "analyzeRecalibratedStarAssociation_gatherResourceUsage_metadata": "json",
              "analyzeRecalibratedStarAssociation_log": "json",
              "analyzeRecalibratedStarAssociation_metadata": "json",
              "analyzeRecalibratedStarObjectMatch_gatherResourceUsage_log": "json",
              "analyzeRecalibratedStarObjectMatch_gatherResourceUsage_metadata": "json",
              "analyzeRecalibratedStarObjectMatch_log": "json",
              "analyzeRecalibratedStarObjectMatch_metadata": "json",
              "analyzeSingleVisitStarAssociation_gatherResourceUsage_log": "json",
              "analyzeSingleVisitStarAssociation_gatherResourceUsage_metadata": "json",
              "analyzeSingleVisitStarAssociation_log": "json",
              "analyzeSingleVisitStarAssociation_metadata": "json",
              "analyzeSourceAssociation_gatherResourceUsage_log": "json",
              "analyzeSourceAssociation_gatherResourceUsage_metadata": "json",
              "analyzeSourceAssociation_log": "json",
              "analyzeSourceAssociation_metadata": "json",
              "assembleDeepCoadd_gatherResourceUsage_log": "json",
              "assembleDeepCoadd_gatherResourceUsage_metadata": "json",
              "assembleDeepCoadd_log": "json",
              "assembleDeepCoadd_metadata": "json",
              "assembleTemplateCoadd_gatherResourceUsage_log": "json",
              "assembleTemplateCoadd_gatherResourceUsage_metadata": "json",
              "assembleTemplateCoadd_log": "json",
              "assembleTemplateCoadd_metadata": "json",
              "associateAnalysisSource_gatherResourceUsage_log": "json",
              "associateAnalysisSource_gatherResourceUsage_metadata": "json",
              "associateAnalysisSource_log": "json",
              "associateAnalysisSource_metadata": "json",
              "associateDiaSource_gatherResourceUsage_log": "json",
              "associateDiaSource_gatherResourceUsage_metadata": "json",
              "associateDiaSource_log": "json",
              "associateDiaSource_metadata": "json",
              "associateIsolatedStar_gatherResourceUsage_log": "json",
              "associateIsolatedStar_gatherResourceUsage_metadata": "json",
              "associateIsolatedStar_log": "json",
              "associateIsolatedStar_metadata": "json",
              "calculateDiaObject_gatherResourceUsage_log": "json",
              "calculateDiaObject_gatherResourceUsage_metadata": "json",
              "calculateDiaObject_log": "json",
              "calculateDiaObject_metadata": "json",
              "calibrateImage_gatherResourceUsage_log": "json",
              "calibrateImage_gatherResourceUsage_metadata": "json",
              "calibrateImage_log": "json",
              "calibrateImage_metadata": "json",
              "catalogMatchTract_gatherResourceUsage_log": "json",
              "catalogMatchTract_gatherResourceUsage_metadata": "json",
              "catalogMatchTract_log": "json",
              "catalogMatchTract_metadata": "json",
              "computeObjectEpochs_gatherResourceUsage_log": "json",
              "computeObjectEpochs_gatherResourceUsage_metadata": "json",
              "computeObjectEpochs_log": "json",
              "computeObjectEpochs_metadata": "json",
              "computeReliability_gatherResourceUsage_log": "json",
              "computeReliability_gatherResourceUsage_metadata": "json",
              "computeReliability_log": "json",
              "computeReliability_metadata": "json",
              "consolidateDiaObject_gatherResourceUsage_log": "json",
              "consolidateDiaObject_gatherResourceUsage_metadata": "json",
              "consolidateDiaObject_log": "json",
              "consolidateDiaObject_metadata": "json",
              "consolidateDiaSource_gatherResourceUsage_log": "json",
              "consolidateDiaSource_gatherResourceUsage_metadata": "json",
              "consolidateDiaSource_log": "json",
              "consolidateDiaSource_metadata": "json",
              "consolidateHealSparsePropertyMaps_gatherResourceUsage_log": "json",
              "consolidateHealSparsePropertyMaps_gatherResourceUsage_metadata": "json",
              "consolidateHealSparsePropertyMaps_log": "json",
              "consolidateHealSparsePropertyMaps_metadata": "json",
              "consolidateObject_gatherResourceUsage_log": "json",
              "consolidateObject_gatherResourceUsage_metadata": "json",
              "consolidateObject_log": "json",
              "consolidateObject_metadata": "json",
              "consolidateRecalibratedStar_gatherResourceUsage_log": "json",
              "consolidateRecalibratedStar_gatherResourceUsage_metadata": "json",
              "consolidateRecalibratedStar_log": "json",
              "consolidateRecalibratedStar_metadata": "json",
              "consolidateResourceUsage_log": "json",
              "consolidateResourceUsage_metadata": "json",
              "consolidateSingleVisitStar_gatherResourceUsage_log": "json",
              "consolidateSingleVisitStar_gatherResourceUsage_metadata": "json",
              "consolidateSingleVisitStar_log": "json",
              "consolidateSingleVisitStar_metadata": "json",
              "consolidateSource_gatherResourceUsage_log": "json",
              "consolidateSource_gatherResourceUsage_metadata": "json",
              "consolidateSource_log": "json",
              "consolidateSource_metadata": "json",
              "consolidateSsTables_gatherResourceUsage_log": "json",
              "consolidateSsTables_gatherResourceUsage_metadata": "json",
              "consolidateSsTables_log": "json",
              "consolidateSsTables_metadata": "json",
              "consolidateVisitDiaSource_gatherResourceUsage_log": "json",
              "consolidateVisitDiaSource_gatherResourceUsage_metadata": "json",
              "consolidateVisitDiaSource_log": "json",
              "consolidateVisitDiaSource_metadata": "json",
              "consolidateVisitSummary_gatherResourceUsage_log": "json",
              "consolidateVisitSummary_gatherResourceUsage_metadata": "json",
              "consolidateVisitSummary_log": "json",
              "consolidateVisitSummary_metadata": "json",
              "deblendCoaddFootprints_gatherResourceUsage_log": "json",
              "deblendCoaddFootprints_gatherResourceUsage_metadata": "json",
              "deblendCoaddFootprints_log": "json",
              "deblendCoaddFootprints_metadata": "json",
              "detectAndMeasureDiaSource_gatherResourceUsage_log": "json",
              "detectAndMeasureDiaSource_gatherResourceUsage_metadata": "json",
              "detectAndMeasureDiaSource_log": "json",
              "detectAndMeasureDiaSource_metadata": "json",
              "detectCoaddPeaks_gatherResourceUsage_log": "json",
              "detectCoaddPeaks_gatherResourceUsage_metadata": "json",
              "detectCoaddPeaks_log": "json",
              "detectCoaddPeaks_metadata": "json",
              "fgcmBuildFromIsolatedStar_gatherResourceUsage_log": "json",
              "fgcmBuildFromIsolatedStar_gatherResourceUsage_metadata": "json",
              "fgcmBuildFromIsolatedStar_log": "json",
              "fgcmBuildFromIsolatedStar_metadata": "json",
              "fgcmFitCycle_gatherResourceUsage_log": "json",
              "fgcmFitCycle_gatherResourceUsage_metadata": "json",
              "fgcmFitCycle_log": "json",
              "fgcmFitCycle_metadata": "json",
              "fgcmOutputProducts_gatherResourceUsage_log": "json",
              "fgcmOutputProducts_gatherResourceUsage_metadata": "json",
              "fgcmOutputProducts_log": "json",
              "fgcmOutputProducts_metadata": "json",
              "filterDiaSource_gatherResourceUsage_log": "json",
              "filterDiaSource_gatherResourceUsage_metadata": "json",
              "filterDiaSource_log": "json",
              "filterDiaSource_metadata": "json",
              "fitDeblendedObjectsSersic_gatherResourceUsage_log": "json",
              "fitDeblendedObjectsSersic_gatherResourceUsage_metadata": "json",
              "fitDeblendedObjectsSersic_log": "json",
              "fitDeblendedObjectsSersic_metadata": "json",
              "fitDeepCoaddPsfGaussians_gatherResourceUsage_log": "json",
              "fitDeepCoaddPsfGaussians_gatherResourceUsage_metadata": "json",
              "fitDeepCoaddPsfGaussians_log": "json",
              "fitDeepCoaddPsfGaussians_metadata": "json",
              "forcedPhotDiaObjectDifference_gatherResourceUsage_log": "json",
              "forcedPhotDiaObjectDifference_gatherResourceUsage_metadata": "json",
              "forcedPhotDiaObjectDifference_log": "json",
              "forcedPhotDiaObjectDifference_metadata": "json",
              "forcedPhotDiaObjectDirect_gatherResourceUsage_log": "json",
              "forcedPhotDiaObjectDirect_gatherResourceUsage_metadata": "json",
              "forcedPhotDiaObjectDirect_log": "json",
              "forcedPhotDiaObjectDirect_metadata": "json",
              "forcedPhotObjectDifference_gatherResourceUsage_log": "json",
              "forcedPhotObjectDifference_gatherResourceUsage_metadata": "json",
              "forcedPhotObjectDifference_log": "json",
              "forcedPhotObjectDifference_metadata": "json",
              "forcedPhotObjectDirect_gatherResourceUsage_log": "json",
              "forcedPhotObjectDirect_gatherResourceUsage_metadata": "json",
              "forcedPhotObjectDirect_log": "json",
              "forcedPhotObjectDirect_metadata": "json",
              "gbdesAstrometricFit_gatherResourceUsage_log": "json",
              "gbdesAstrometricFit_gatherResourceUsage_metadata": "json",
              "gbdesAstrometricFit_log": "json",
              "gbdesAstrometricFit_metadata": "json",
              "isr_gatherResourceUsage_log": "json",
              "isr_gatherResourceUsage_metadata": "json",
              "isr_log": "json",
              "isr_metadata": "json",
              "makeAnalysisRecalibratedStarAssociationMetricTable_gatherResourceUsage_log": "json",
              "makeAnalysisRecalibratedStarAssociationMetricTable_gatherResourceUsage_metadata": "json",
              "makeAnalysisRecalibratedStarAssociationMetricTable_log": "json",
              "makeAnalysisRecalibratedStarAssociationMetricTable_metadata": "json",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_gatherResourceUsage_log": "json",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_gatherResourceUsage_metadata": "json",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_log": "json",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_metadata": "json",
              "makeAnalysisSingleVisitStarAssociationMetricTable_gatherResourceUsage_log": "json",
              "makeAnalysisSingleVisitStarAssociationMetricTable_gatherResourceUsage_metadata": "json",
              "makeAnalysisSingleVisitStarAssociationMetricTable_log": "json",
              "makeAnalysisSingleVisitStarAssociationMetricTable_metadata": "json",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_gatherResourceUsage_log": "json",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_gatherResourceUsage_metadata": "json",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_log": "json",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_metadata": "json",
              "makeAnalysisSourceAssociationMetricTable_gatherResourceUsage_log": "json",
              "makeAnalysisSourceAssociationMetricTable_gatherResourceUsage_metadata": "json",
              "makeAnalysisSourceAssociationMetricTable_log": "json",
              "makeAnalysisSourceAssociationMetricTable_metadata": "json",
              "makeAnalysisSourceAssociationWholeSkyPlot_gatherResourceUsage_config": "json",
              "makeAnalysisSourceAssociationWholeSkyPlot_gatherResourceUsage_log": "json",
              "makeAnalysisSourceAssociationWholeSkyPlot_gatherResourceUsage_metadata": "json",
              "makeAnalysisSourceAssociationWholeSkyPlot_log": "json",
              "makeAnalysisSourceAssociationWholeSkyPlot_metadata": "json",
              "makeDirectWarp_gatherResourceUsage_log": "json",
              "makeDirectWarp_gatherResourceUsage_metadata": "json",
              "makeDirectWarp_log": "json",
              "makeDirectWarp_metadata": "json",
              "makeHealSparsePropertyMaps_gatherResourceUsage_log": "json",
              "makeHealSparsePropertyMaps_gatherResourceUsage_metadata": "json",
              "makeHealSparsePropertyMaps_log": "json",
              "makeHealSparsePropertyMaps_metadata": "json",
              "makeInitialVisitDetectorTable_gatherResourceUsage_log": "json",
              "makeInitialVisitDetectorTable_gatherResourceUsage_metadata": "json",
              "makeInitialVisitDetectorTable_log": "json",
              "makeInitialVisitDetectorTable_metadata": "json",
              "makeInitialVisitTable_gatherResourceUsage_log": "json",
              "makeInitialVisitTable_gatherResourceUsage_metadata": "json",
              "makeInitialVisitTable_log": "json",
              "makeInitialVisitTable_metadata": "json",
              "makeMetricTableObjectTableCoreRefCatMatch_gatherResourceUsage_log": "json",
              "makeMetricTableObjectTableCoreRefCatMatch_gatherResourceUsage_metadata": "json",
              "makeMetricTableObjectTableCoreRefCatMatch_log": "json",
              "makeMetricTableObjectTableCoreRefCatMatch_metadata": "json",
              "makeMetricTableObjectTableCore_gatherResourceUsage_log": "json",
              "makeMetricTableObjectTableCore_gatherResourceUsage_metadata": "json",
              "makeMetricTableObjectTableCore_log": "json",
              "makeMetricTableObjectTableCore_metadata": "json",
              "makePsfMatchedWarp_gatherResourceUsage_log": "json",
              "makePsfMatchedWarp_gatherResourceUsage_metadata": "json",
              "makePsfMatchedWarp_log": "json",
              "makePsfMatchedWarp_metadata": "json",
              "makeVisitDetectorTable_gatherResourceUsage_log": "json",
              "makeVisitDetectorTable_gatherResourceUsage_metadata": "json",
              "makeVisitDetectorTable_log": "json",
              "makeVisitDetectorTable_metadata": "json",
              "makeVisitTable_gatherResourceUsage_log": "json",
              "makeVisitTable_gatherResourceUsage_metadata": "json",
              "makeVisitTable_log": "json",
              "makeVisitTable_metadata": "json",
              "measureObjectForced_gatherResourceUsage_log": "json",
              "measureObjectForced_gatherResourceUsage_metadata": "json",
              "measureObjectForced_log": "json",
              "measureObjectForced_metadata": "json",
              "measureObjectUnforced_gatherResourceUsage_log": "json",
              "measureObjectUnforced_gatherResourceUsage_metadata": "json",
              "measureObjectUnforced_log": "json",
              "measureObjectUnforced_metadata": "json",
              "mergeObjectDetection_gatherResourceUsage_log": "json",
              "mergeObjectDetection_gatherResourceUsage_metadata": "json",
              "mergeObjectDetection_log": "json",
              "mergeObjectDetection_metadata": "json",
              "mergeObjectMeasurement_gatherResourceUsage_log": "json",
              "mergeObjectMeasurement_gatherResourceUsage_metadata": "json",
              "mergeObjectMeasurement_log": "json",
              "objectTableCoreRefCatMatchWholeSkyPlot_gatherResourceUsage_log": "json",
              "objectTableCoreRefCatMatchWholeSkyPlot_gatherResourceUsage_metadata": "json",
              "objectTableCoreRefCatMatchWholeSkyPlot_log": "json",
              "objectTableCoreRefCatMatchWholeSkyPlot_metadata": "json",
              "objectTableCoreWholeSkyPlot_gatherResourceUsage_log": "json",
              "objectTableCoreWholeSkyPlot_gatherResourceUsage_metadata": "json",
              "objectTableCoreWholeSkyPlot_log": "json",
              "objectTableCoreWholeSkyPlot_metadata": "json",
              "object_scarlet_models": "json",
              "photometricCatalogMatch_gatherResourceUsage_log": "json",
              "photometricCatalogMatch_gatherResourceUsage_metadata": "json",
              "photometricCatalogMatch_log": "json",
              "photometricCatalogMatch_metadata": "json",
              "photometricRefCatObjectTract_gatherResourceUsage_log": "json",
              "photometricRefCatObjectTract_gatherResourceUsage_metadata": "json",
              "photometricRefCatObjectTract_log": "json",
              "photometricRefCatObjectTract_metadata": "json",
              "plotPropertyMapSurvey_gatherResourceUsage_log": "json",
              "plotPropertyMapSurvey_gatherResourceUsage_metadata": "json",
              "plotPropertyMapSurvey_log": "json",
              "plotPropertyMapSurvey_metadata": "json",
              "plotPropertyMapTract_gatherResourceUsage_log": "json",
              "plotPropertyMapTract_gatherResourceUsage_metadata": "json",
              "plotPropertyMapTract_log": "json",
              "plotPropertyMapTract_metadata": "json",
              "recalibrateSingleVisitStar_gatherResourceUsage_log": "json",
              "recalibrateSingleVisitStar_gatherResourceUsage_metadata": "json",
              "recalibrateSingleVisitStar_log": "json",
              "recalibrateSingleVisitStar_metadata": "json",
              "refCatObjectTract_gatherResourceUsage_log": "json",
              "refCatObjectTract_gatherResourceUsage_metadata": "json",
              "refCatObjectTract_log": "json",
              "refitPsfModel_gatherResourceUsage_log": "json",
              "refitPsfModel_gatherResourceUsage_metadata": "json",
              "refitPsfModel_metadata": "json",
              "reprocessVisitImage_gatherResourceUsage_log": "json",
              "reprocessVisitImage_gatherResourceUsage_metadata": "json",
              "reprocessVisitImage_log": "json",
              "reprocessVisitImage_metadata": "json",
              "rewarpTemplate_gatherResourceUsage_log": "json",
              "rewarpTemplate_gatherResourceUsage_metadata": "json",
              "rewarpTemplate_log": "json",
              "rewarpTemplate_metadata": "json",
              "rewriteObject_gatherResourceUsage_log": "json",
              "rewriteObject_gatherResourceUsage_metadata": "json",
              "rewriteObject_log": "json",
              "rewriteObject_metadata": "json",
              "selectDeepCoaddVisits_gatherResourceUsage_log": "json",
              "selectDeepCoaddVisits_gatherResourceUsage_metadata": "json",
              "selectDeepCoaddVisits_log": "json",
              "selectDeepCoaddVisits_metadata": "json",
              "selectTemplateCoaddVisits_gatherResourceUsage_log": "json",
              "selectTemplateCoaddVisits_gatherResourceUsage_metadata": "json",
              "selectTemplateCoaddVisits_log": "json",
              "selectTemplateCoaddVisits_metadata": "json",
              "splitPrimaryObjectForcedSource_gatherResourceUsage_log": "json",
              "splitPrimaryObjectForcedSource_gatherResourceUsage_metadata": "json",
              "splitPrimaryObjectForcedSource_log": "json",
              "splitPrimaryObjectForcedSource_metadata": "json",
              "splitPrimaryObject_gatherResourceUsage_log": "json",
              "splitPrimaryObject_gatherResourceUsage_metadata": "json",
              "splitPrimaryObject_log": "json",
              "splitPrimaryObject_metadata": "json",
              "splitPrimarySource_gatherResourceUsage_log": "json",
              "splitPrimarySource_gatherResourceUsage_metadata": "json",
              "splitPrimarySource_log": "json",
              "splitPrimarySource_metadata": "json",
              "standardizeDiaObjectForcedSource_gatherResourceUsage_log": "json",
              "standardizeDiaObjectForcedSource_gatherResourceUsage_metadata": "json",
              "standardizeDiaObjectForcedSource_log": "json",
              "standardizeDiaObjectForcedSource_metadata": "json",
              "standardizeDiaSource_gatherResourceUsage_log": "json",
              "standardizeDiaSource_gatherResourceUsage_metadata": "json",
              "standardizeDiaSource_log": "json",
              "standardizeDiaSource_metadata": "json",
              "standardizeObjectForcedSource_gatherResourceUsage_log": "json",
              "standardizeObjectForcedSource_gatherResourceUsage_metadata": "json",
              "standardizeObjectForcedSource_log": "json",
              "standardizeObjectForcedSource_metadata": "json",
              "standardizeObject_gatherResourceUsage_log": "json",
              "standardizeObject_gatherResourceUsage_metadata": "json",
              "standardizeObject_log": "json",
              "standardizeObject_metadata": "json",
              "standardizeRecalibratedStar_gatherResourceUsage_log": "json",
              "standardizeRecalibratedStar_gatherResourceUsage_metadata": "json",
              "standardizeRecalibratedStar_log": "json",
              "standardizeRecalibratedStar_metadata": "json",
              "standardizeSingleVisitStar_gatherResourceUsage_log": "json",
              "standardizeSingleVisitStar_gatherResourceUsage_metadata": "json",
              "standardizeSingleVisitStar_log": "json",
              "standardizeSingleVisitStar_metadata": "json",
              "standardizeSource_gatherResourceUsage_log": "json",
              "standardizeSource_gatherResourceUsage_metadata": "json",
              "standardizeSource_log": "json",
              "standardizeSource_metadata": "json",
              "subtractImages_gatherResourceUsage_log": "json",
              "subtractImages_gatherResourceUsage_metadata": "json",
              "subtractImages_log": "json",
              "subtractImages_metadata": "json",
              "updateVisitSummary_gatherResourceUsage_log": "json",
              "updateVisitSummary_gatherResourceUsage_metadata": "json",
              "updateVisitSummary_log": "json",
              "updateVisitSummary_metadata": "json",
              "validateObjectTableCore_gatherResourceUsage_log": "json",
              "validateObjectTableCore_gatherResourceUsage_metadata": "json",
              "validateObjectTableCore_log": "json",
              "validateObjectTableCore_metadata": "json",
              "writeDiaObjectForcedSource_gatherResourceUsage_log": "json",
              "writeDiaObjectForcedSource_gatherResourceUsage_metadata": "json",
              "writeDiaObjectForcedSource_log": "json",
              "writeDiaObjectForcedSource_metadata": "json",
              "writeObjectForcedSource_gatherResourceUsage_log": "json",
              "writeObjectForcedSource_gatherResourceUsage_metadata": "json",
              "writeObjectForcedSource_log": "json",
              "writeObjectForcedSource_metadata": "json",
              "refitPsfModel_log": "json",         
              "makeInitialVisitDetectorTable_config": "py",
              "analyzeDiaSourceTableTract_config": "py",
              "analyzeDiaSourceTableTract_gatherResourceUsage_config": "py",
              "analyzeObjectTableCore_config": "py",
              "analyzeObjectTableCore_gatherResourceUsage_config": "py",
              "analyzeObjectTableSurveyCore_config": "py",
              "analyzeObjectTableSurveyCore_gatherResourceUsage_config": "py",
              "analyzeRecalibratedStarAssociation_config": "py",
              "analyzeRecalibratedStarAssociation_gatherResourceUsage_config": "py",
              "analyzeRecalibratedStarObjectMatch_config": "py",
              "analyzeRecalibratedStarObjectMatch_gatherResourceUsage_config": "py",
              "analyzeSingleVisitStarAssociation_config": "py",
              "analyzeSingleVisitStarAssociation_gatherResourceUsage_config": "py",
              "analyzeSourceAssociation_config": "py",
              "analyzeSourceAssociation_gatherResourceUsage_config": "py",
              "assembleDeepCoadd_config": "py",
              "assembleDeepCoadd_gatherResourceUsage_config": "py",
              "assembleTemplateCoadd_config": "py",
              "assembleTemplateCoadd_gatherResourceUsage_config": "py",
              "associateAnalysisSource_config": "py",
              "associateAnalysisSource_gatherResourceUsage_config": "py",
              "associateDiaSource_config": "py",
              "associateDiaSource_gatherResourceUsage_config": "py",
              "associateIsolatedStar_config": "py",
              "associateIsolatedStar_gatherResourceUsage_config": "py",
              "calculateDiaObject_config": "py",
              "calculateDiaObject_gatherResourceUsage_config": "py",
              "calibrateImage_config": "py",
              "calibrateImage_gatherResourceUsage_config": "py",
              "catalogMatchTract_config": "py",
              "catalogMatchTract_gatherResourceUsage_config": "py",
              "computeObjectEpochs_config": "py",
              "computeObjectEpochs_gatherResourceUsage_config": "py",
              "computeReliability_config": "py",
              "computeReliability_gatherResourceUsage_config": "py",
              "consolidateDiaObject_config": "py",
              "consolidateDiaObject_gatherResourceUsage_config": "py",
              "consolidateDiaSource_config": "py",
              "consolidateDiaSource_gatherResourceUsage_config": "py",
              "consolidateHealSparsePropertyMaps_config": "py",
              "consolidateHealSparsePropertyMaps_gatherResourceUsage_config": "py",
              "consolidateObject_config": "py",
              "consolidateObject_gatherResourceUsage_config": "py",
              "consolidateRecalibratedStar_config": "py",
              "consolidateRecalibratedStar_gatherResourceUsage_config": "py",
              "consolidateResourceUsage_config": "py",
              "consolidateSingleVisitStar_config": "py",
              "consolidateSingleVisitStar_gatherResourceUsage_config": "py",
              "consolidateSource_config": "py",
              "consolidateSource_gatherResourceUsage_config": "py",
              "consolidateSsTables_config": "py",
              "consolidateSsTables_gatherResourceUsage_config": "py",
              "consolidateVisitDiaSource_config": "py",
              "consolidateVisitDiaSource_gatherResourceUsage_config": "py",
              "consolidateVisitSummary_config": "py",
              "consolidateVisitSummary_gatherResourceUsage_config": "py",
              "deblendCoaddFootprints_config": "py",
              "deblendCoaddFootprints_gatherResourceUsage_config": "py",
              "detectAndMeasureDiaSource_config": "py",
              "detectAndMeasureDiaSource_gatherResourceUsage_config": "py",
              "detectCoaddPeaks_config": "py",
              "detectCoaddPeaks_gatherResourceUsage_config": "py",
              "fgcmBuildFromIsolatedStar_config": "py",
              "fgcmBuildFromIsolatedStar_gatherResourceUsage_config": "py",
              "fgcmFitCycle_config": "py",
              "fgcmFitCycle_gatherResourceUsage_config": "py",
              "fgcmOutputProducts_config": "py",
              "fgcmOutputProducts_gatherResourceUsage_config": "py",
              "filterDiaSource_config": "py",
              "filterDiaSource_gatherResourceUsage_config": "py",
              "fitDeblendedObjectsSersic_config": "py",
              "fitDeblendedObjectsSersic_gatherResourceUsage_config": "py",
              "fitDeepCoaddPsfGaussians_config": "py",
              "fitDeepCoaddPsfGaussians_gatherResourceUsage_config": "py",
              "forcedPhotDiaObjectDifference_config": "py",
              "forcedPhotDiaObjectDifference_gatherResourceUsage_config": "py",
              "forcedPhotDiaObjectDirect_config": "py",
              "forcedPhotDiaObjectDirect_gatherResourceUsage_config": "py",
              "forcedPhotObjectDifference_config": "py",
              "forcedPhotObjectDifference_gatherResourceUsage_config": "py",
              "forcedPhotObjectDirect_config": "py",
              "forcedPhotObjectDirect_gatherResourceUsage_config": "py",
              "gbdesAstrometricFit_config": "py",
              "gbdesAstrometricFit_gatherResourceUsage_config": "py",
              "isr_config": "py",
              "isr_gatherResourceUsage_config": "py",
              "makeAnalysisRecalibratedStarAssociationMetricTable_config": "py",
              "makeAnalysisRecalibratedStarAssociationMetricTable_gatherResourceUsage_config": "py",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_config": "py",
              "makeAnalysisRecalibratedStarAssociationWholeSkyPlot_gatherResourceUsage_config": "py",
              "makeAnalysisSingleVisitStarAssociationMetricTable_config": "py",
              "makeAnalysisSingleVisitStarAssociationMetricTable_gatherResourceUsage_config": "py",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_config": "py",
              "makeAnalysisSingleVisitStarAssociationWholeSkyPlot_gatherResourceUsage_config": "py",
              "makeAnalysisSourceAssociationMetricTable_config": "py",
              "makeAnalysisSourceAssociationMetricTable_gatherResourceUsage_config": "py",
              "makeAnalysisSourceAssociationWholeSkyPlot_config": "py",
              "makeDirectWarp_config": "py",
              "makeDirectWarp_gatherResourceUsage_config": "py",
              "makeHealSparsePropertyMaps_config": "py",
              "makeHealSparsePropertyMaps_gatherResourceUsage_config": "py",
              "makeInitialVisitDetectorTable_gatherResourceUsage_config": "py",
              "makeInitialVisitTable_config": "py",
              "makeInitialVisitTable_gatherResourceUsage_config": "py",
              "makeMetricTableObjectTableCoreRefCatMatch_config": "py",
              "makeMetricTableObjectTableCoreRefCatMatch_gatherResourceUsage_config": "py",
              "makeMetricTableObjectTableCore_config": "py",
              "makeMetricTableObjectTableCore_gatherResourceUsage_config": "py",
              "makePsfMatchedWarp_config": "py",
              "makePsfMatchedWarp_gatherResourceUsage_config": "py",
              "makeVisitDetectorTable_config": "py",
              "makeVisitDetectorTable_gatherResourceUsage_config": "py",
              "makeVisitTable_config": "py",
              "makeVisitTable_gatherResourceUsage_config": "py",
              "measureObjectForced_config": "py",
              "measureObjectForced_gatherResourceUsage_config": "py",
              "measureObjectUnforced_config": "py",
              "measureObjectUnforced_gatherResourceUsage_config": "py",
              "mergeObjectDetection_config": "py",
              "mergeObjectDetection_gatherResourceUsage_config": "py",
              "mergeObjectMeasurement_config": "py",
              "mergeObjectMeasurement_gatherResourceUsage_config": "py",
              "objectTableCoreRefCatMatchWholeSkyPlot_config": "py",
              "objectTableCoreRefCatMatchWholeSkyPlot_gatherResourceUsage_config": "py",
              "objectTableCoreWholeSkyPlot_config": "py",
              "objectTableCoreWholeSkyPlot_gatherResourceUsage_config": "py",
              "photometricCatalogMatch_config": "py",
              "photometricCatalogMatch_gatherResourceUsage_config": "py",
              "photometricRefCatObjectTract_config": "py",
              "photometricRefCatObjectTract_gatherResourceUsage_config": "py",
              "plotPropertyMapSurvey_config": "py",
              "plotPropertyMapSurvey_gatherResourceUsage_config": "py",
              "plotPropertyMapTract_config": "py",
              "plotPropertyMapTract_gatherResourceUsage_config": "py",
              "recalibrateSingleVisitStar_config": "py",
              "recalibrateSingleVisitStar_gatherResourceUsage_config": "py",
              "refCatObjectTract_config": "py",
              "refCatObjectTract_gatherResourceUsage_config": "py",
              "refitPsfModel_config": "py",
              "refitPsfModel_gatherResourceUsage_config": "py",
              "reprocessVisitImage_config": "py",
              "reprocessVisitImage_gatherResourceUsage_config": "py",
              "rewarpTemplate_config": "py",
              "rewarpTemplate_gatherResourceUsage_config": "py",
              "rewriteObject_config": "py",
              "rewriteObject_gatherResourceUsage_config": "py",
              "selectDeepCoaddVisits_config": "py",
              "selectDeepCoaddVisits_gatherResourceUsage_config": "py",
              "selectTemplateCoaddVisits_config": "py",
              "selectTemplateCoaddVisits_gatherResourceUsage_config": "py",
              "splitPrimaryObjectForcedSource_config": "py",
              "splitPrimaryObjectForcedSource_gatherResourceUsage_config": "py",
              "splitPrimaryObject_config": "py",
              "splitPrimaryObject_gatherResourceUsage_config": "py",
              "splitPrimarySource_config": "py",
              "splitPrimarySource_gatherResourceUsage_config": "py",
              "standardizeDiaObjectForcedSource_config": "py",
              "standardizeDiaObjectForcedSource_gatherResourceUsage_config": "py",
              "standardizeDiaSource_config": "py",
              "standardizeDiaSource_gatherResourceUsage_config": "py",
              "standardizeObjectForcedSource_config": "py",
              "standardizeObjectForcedSource_gatherResourceUsage_config": "py",
              "standardizeObject_config": "py",
              "standardizeObject_gatherResourceUsage_config": "py",
              "standardizeRecalibratedStar_config": "py",
              "standardizeRecalibratedStar_gatherResourceUsage_config": "py",
              "standardizeSingleVisitStar_config": "py",
              "standardizeSingleVisitStar_gatherResourceUsage_config": "py",
              "standardizeSource_config": "py",
              "standardizeSource_gatherResourceUsage_config": "py",
              "subtractImages_config": "py",
              "subtractImages_gatherResourceUsage_config": "py",
              "updateVisitSummary_config": "py",
              "updateVisitSummary_gatherResourceUsage_config": "py",
              "validateObjectTableCore_config": "py",
              "validateObjectTableCore_gatherResourceUsage_config": "py",
              "writeDiaObjectForcedSource_config": "py",
              "writeDiaObjectForcedSource_gatherResourceUsage_config": "py",
              "writeObjectForcedSource_config": "py",
              "writeObjectForcedSource_gatherResourceUsage_config": "py"
            }
detectorIdComCam = {"R22_S00": 0, "R22_S01": 1, "R22_S02": 2,
                    "R22_S10": 3, "R22_S11": 4, "R22_S12": 5,
                    "R22_S20": 6, "R22_S21": 7, "R22_S22": 8}
verbose = 0

# set panda options for table views
pd.options.display.max_rows = None
pd.options.display.max_columns = None
pd.options.display.max_colwidth = None
pd.options.display.expand_frame_repr = False

#-----------------------------------------------------------------------------

def change_base_directory(url):
    directories=url.split('/')[5:]
    filename = url.rsplit('/',1)[1]

    if "raw" in directories[0]:
        dtype = directories[2]
    elif "DM" in directories[5]:
        dtype = directories[6]
    elif "v29" in directories[5]:
        dtype = directories[8]
    if verbose > 10:
        print(directories,":::",dtype)

    if dtype == "deep_coadd":
        band = directories[-2]
        patch = int(directories[-3])
        tract = int(directories[-4])
        run = directories[-6]
    elif dtype == "template_coadd":
        band = directories[-2]
        patch = int(directories[-3])
        tract = int(directories[-4])
        run = directories[-6]
    elif dtype == "visit_image":
        band = directories[-4]
        instrument = "LSSTComCam"
        day_obs = int(directories[-5])
        detector = detectorIdComCam['_'.join(filename.split('_')[7:9])]
        physical_filter =directories[-3]
        visit = int(directories[-2])
    elif dtype == "difference_image":
        band = directories[-4]
        instrument = "LSSTComCam"
        day_obs = int(directories[-5])
        detector = detectorIdComCam['_'.join(filename.split('_')[7:9])]
        physical_filter =directories[-3]
        visit = int(directories[-2])
    elif dtype == "object":
        tract = int(directories[-2])
    elif dtype == "source":
        band = directories[-4]
        instrument = "LSSTComCam"
        day_obs = int(directories[-5])
        physical_filter = directories[-3]
        visit = int(directories[-2])
    elif dtype == "object_forced_source":
        tract = int(directories[-3])
        patch = int(directories[-2])
    elif dtype == "dia_object":
        tract = int(directories[-2])
    elif dtype == "dia_source":
        tract = int(directories[-2])
    elif dtype == "dia_object_forced_source":
        tract = int(directories[-3])
        patch = int(directories[-2])
    elif dtype == "ss_object":
        instrument = "LSSTComCam"
    elif dtype == "ss_source":
        instrument = "LSSTComCam"
    elif dtype == "visit_table":
        instrument = "LSSTComCam"
    elif dtype == "visit_detector_table":
        instrument = "LSSTComCam"
    elif dtype == "deep_coadd_background":
        band = directories[-2]
        tract = int(directories[-4])
        patch = int(directories[-3])
    elif dtype == "deep_coadd_n_image":
        band = directories[-2]
        tract = int(directories[-4])
        patch = int(directories[-3])
    elif dtype == "visit_summary":
        band = directories[-4]
        instrument = "LSSTComCam"
        day_obs = int(directories[-5])
        physical_filter = directories[-3]
        visit = int(directories[-2])
    elif dtype == "visit_image_background":
        band = directories[-4]
        instrument = "LSSTComCam"
        day_obs = int(directories[-5])
        detector = detectorIdComCam['_'.join(filename.split('_')[8:10])]
        physical_filter = directories[-3]
        visit = int(directories[-2])
    elif dtype == "the_monster_20250219":
        htm7 = int(filename.split('.')[0])
    elif dtype == "deepCoadd_dcr_ddec_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_dcr_dra_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_dcr_e1_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_dcr_e2_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_epoch_consolidated_map_max":
        band = directories[-2]
    elif dtype == "deepCoadd_epoch_consolidated_map_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_epoch_consolidated_map_min":
        band = directories[-2]
    elif dtype == "deepCoadd_exposure_time_consolidated_map_sum":
        band = directories[-2]
    elif dtype == "deepCoadd_psf_e1_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_psf_e2_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_psf_maglim_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_psf_size_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_sky_background_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "deepCoadd_sky_noise_consolidated_map_weighted_mean":
        band = directories[-2]
    elif dtype == "mergeObjectMeasurement_metadata":
        tract = int(directories[-3])
        patch = int(directories[-2])
    elif dtype == "refCatObjectTract_metadata":
        tract = int(directories[-2])
    else:
        print(f"[ERROR] dataType {dtype} not supported.")
        raise SystemExit

    skymap = 'lsst_cells_v1'

    new_url = os.path.join(baseDir, *directories)
    if verbose > 10:
        print("URL:",url)
        print("DIR:",len(directories),directories)
        print("NURL:",new_url)
        #print("bstpr:",band, skymap, tract, patch, run)
        raise SystemExit
    if "coadd" in dtype:
        return new_url, band, skymap, tract, patch, run
    elif "image" in dtype:
        return new_url, band, instrument, day_obs, detector, physical_filter, visit
    elif dtype == "object":
        return new_url, skymap, tract
    elif dtype == "source":
        return new_url, band, instrument, day_obs, physical_filter, visit
    elif dtype == "object_forced_source":
        return new_url, skymap, tract, patch
    elif dtype == "dia_object":
        return new_url, skymap, tract
    elif dtype == "dia_source":
        return new_url, skymap, tract
    elif dtype == "dia_object_forced_source":
        return new_url, skymap, tract, patch
    elif dtype == "ss_object":
        return new_url, instrument
    elif dtype == "ss_source":
        return new_url, instrument
    elif dtype == "visit_table":
        return new_url, instrument
    elif dtype == "visit_detector_table":
        return new_url, instrument
    elif dtype == "deep_coadd_background":
        return new_url, band, tract, patch
    elif dtype == "deep_coadd_n_image":
        return new_url, band, tract, patch
    elif dtype == "visit_summary":
        return new_url, band, instrument, day_obs, physical_filter, visit
    elif dtype == "visit_image_background":
        return new_url, band, instrument, day_obs, detector, physical_filter, visit
    elif dtype == "the_monster_20250219":
        return new_url, htm7
    elif dtype == "deepCoadd_dcr_ddec_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_dcr_dra_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_dcr_e1_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_dcr_e2_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_epoch_consolidated_map_max":
        return new_url, band, skymap
    elif dtype == "deepCoadd_epoch_consolidated_map_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_epoch_consolidated_map_min":
        return new_url, band, skymap
    elif dtype == "deepCoadd_exposure_time_consolidated_map_sum":
        return new_url, band, skymap
    elif dtype == "deepCoadd_psf_e1_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_psf_e2_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_psf_maglim_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_psf_size_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_sky_background_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "deepCoadd_sky_noise_consolidated_map_weighted_mean":
        return new_url, band, skymap
    elif dtype == "mergeObjectMeasurement_metadata":
        return new_url, skymap, tract, patch
    elif dtype == "refCatObjectTract_metadata":
        return new_url, skymap, tract
#-----------------------------------------------------------------------------

def read_dirs(dataDir, dataType, fileType):
    csvPath = os.path.join(csvBaseDir, f'{dataType}_urls.csv')
    #globPath = os.path.join(base_dir, 'runs/**/*.fits')
    #globPath = os.path.join(base_dir, 'runs/**/*.parq')
    #globPath = os.path.join(base_dir, 'runs/**/*.%s' % fileType)
    globPath = os.path.join(baseDir, dataDir, '**/*.%s' % fileType)
    print(f"{dataType}: Globbing {globPath} ...")

    filepaths = glob.glob(globPath, recursive=True)
    print("Writing csv file...")
    with open(csvPath, 'w') as f:
        for path in filepaths:
            if dataType == "other":
                exclKeys = list(groups.keys()) + list(excludedDT.keys())
                if all(['/%s/' % dt not in path for dt in exclKeys]):
                    thedtype = path.split('/')[13]
                    theftype = path.rsplit('.',1)[1]
                    f.write(f"              \"{thedtype}\": \"{theftype}\",\n")
                    #f.write(f"{path}\n")
            elif '/%s/' % dataType in path:
                f.write(f"{path}\n")
    print(f"File list written to: {csvPath}")
    
#-----------------------------------------------------------------------------

def processData(dataType, makeList=False, dataDir=defDataDir, fileType='*',
                splitObs=False, removeDupl=False, testRun=False):

    # crawl data dirs and create input url list
    if makeList:
        read_dirs(dataDir, dataType, fileType)
        raise SystemExit

    # read dataType files into 'urls'
    csvPath = os.path.join(csvBaseDir, f'{dataType}_urls.csv')
    urls = pd.read_csv(csvPath, header=None, names=['urls'])

    # change the html link into the base dir
    if "coadd" in dataType:
        urls['filename'], urls['band'], urls['skymap'], urls['tract'], urls['patch'], urls['run'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif "image" in dataType:
        urls['filename'], urls['band'], urls['instrument'], urls['day_obs'], urls['detector'], urls['physical_filter'], urls['visit'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "object":
        urls['filename'], urls['skymap'], urls['tract'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "source":
        urls['filename'], urls['band'], urls['instrument'], urls['day_obs'], urls['physical_filter'], urls['visit'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "object_forced_source":
        urls['filename'], urls['skymap'], urls['tract'], urls['patch'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "dia_object":
        urls['filename'], urls['skymap'], urls['tract'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "dia_source":
        urls['filename'], urls['skymap'], urls['tract'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "dia_object_forced_source":
        urls['filename'], urls['skymap'], urls['tract'], urls['patch'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "ss_object":
        print(urls, urls.shape)
        print(*urls['urls'].apply(lambda x: change_base_directory(x)))
        print(list(zip(*urls['urls'].apply(lambda x: change_base_directory(x)))))
        urls['filename'], urls['instrument'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "ss_source":
        urls['filename'], urls['instrument'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "visit_table":
        urls['filename'], urls['instrument'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "visit_detector_table":
        urls['filename'], urls['instrument'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deep_coadd_background":
        urls['filename'], urls['band'], urls['tract'], urls['patch'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deep_coadd_n_image":
        urls['filename'], urls['band'], urls['tract'], urls['patch'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "visit_summary":
        urls['filename'], urls['band'], urls['instrument'], urls['day_obs'], urls['physical_filter'], urls['visit'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "visit_image_background":
        urls['filename'], urls['band'], urls['instrument'], urls['day_obs'], urls['detector'], urls['physical_filter'], urls['visit'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "the_monster_20250219":
        urls['filename'], urls['htm7'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_dcr_ddec_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_dcr_dra_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_dcr_e1_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_dcr_e2_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_epoch_consolidated_map_max":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_epoch_consolidated_map_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_epoch_consolidated_map_min":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_exposure_time_consolidated_map_sum":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_psf_e1_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_psf_e2_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_psf_maglim_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_psf_size_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_sky_background_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "deepCoadd_sky_noise_consolidated_map_weighted_mean":
        urls['filename'], urls['band'], urls['skymap'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "mergeObjectMeasurement_metadata":
        urls['filename'], urls['skymap'], urls['tract'], urls['patch'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))
    elif dataType == "refCatObjectTract_metadata":
        urls['filename'], urls['skymap'], urls['tract'] = zip(
            *urls['urls'].apply(lambda x: change_base_directory(x)))

        
    if removeDupl:
        # deep copy the existing data to manipulate and delete duplicates
        urlstoo = urls.copy(deep=True)

        # create a duplicated column
        if "coadd" in dataType:
            urlstoo['dupl'] = urlstoo.duplicated(
                subset=['band','tract','patch'], keep='last')
            # sort the duplicates
            #urlstoo.sort_values(by=['band','tract','patch','run'],axis=0,inplace=True)
        elif "image" in dataType:
            urlstoo['dupl'] = urlstoo.duplicated(
                subset=['instrument', 'detector', 'visit'], keep='last')

        if verbose > 9:
            # check the 'urlstoo' object
            print(urlstoo.head())
            if "coadd" in dataType:
                print(urlstoo[['run','band','tract','patch','dupl','filename']].loc[urlstoo['dupl'] == True])
                print(urls['filename'],urls['band'],urls['tract'],urls['patch'])
                print(urls['run'].value_counts())
            elif "image" in dataType:
                print(urlstoo[['instrument','detector','visit','band','day_obs','physical_filter','dupl','filename']].loc[urlstoo['dupl'] == True])
                print(urls['filename'],urls['instrument'],urls['detector'],urls['visit'])
                print(urls['visit'].value_counts())
            raise SystemExit
        # drop the first duplicates and keep the last
        print("w/dupl:", urlstoo.shape)
        if "coadd" in dataType:
            urlstoo.drop_duplicates(subset=['band','tract','patch'],
                                    keep='last', inplace=True,
                                    ignore_index=True)
        elif "image" in dataType:
            urlstoo.drop_duplicates(subset=['instrument', 'detector', 'visit'],
                                    keep='last', inplace=True,
                                    ignore_index=True)
        print("uniqued:", urlstoo.shape)

        # copy the urls list to be used in the file creation
        useUrls = urlstoo.copy(deep=True)
    else:
        useUrls = urls.copy(deep=True)


    try:
        outGroups = groups[dataType]
    except KeyError:
        print("[ERROR] dataType doesn't exist in groups.")
        raise SystemExit

    # create files split by observation run
    if splitObs:
        by_run = useUrls.groupby('run')
        for run_value, group_df in by_run:
            new_group = group_df[outGroups]
            file_table = Table.from_pandas(new_group)
            if verbose > 2:
                print(file_table)

            # set out put path for ecsv file
            ecsvPath = os.path.join(ecsvBaseDir, f'{dataType}_{run_value}.ecsv')

            if testRun:
                print(f"[TEST] write to obsrun tables: {ecsvPath}")
            else:
                print(f"Writing to obsrun table: {ecsvPath}")
                file_table.write(ecsvPath)

    # create one file containing all data
    new_urls = useUrls[outGroups]
    file_table = Table.from_pandas(new_urls)
    if verbose > 2:
        print(file_table)

    # set output path for ecsv file
    ecsvPath = os.path.join(ecsvBaseDir, f'{dataType}.ecsv')

    if testRun:
        print(f"[TEST] write to single table: {ecsvPath}")
    else:
        print(f"Writing to single table: {ecsvPath}")
        file_table.write(ecsvPath)

        
#-----------------------------------------------------------------------------

def usage():
        print("Usage: python parse_to_ecsv.py [-s/--splitobs] [-h/--help] [-t/--test] [-f/--filetype <fileext>] [-l/--makelist] datasetType")
        print("-f/--filetype <fileext>: search for this file extension")
        print("-l/--makelist: create ingest files list")
        print("-d/datadir: top level directory for the ingest files")
        print("-h/--help: print the help")
        print("-r/--remdupl: remove duplicates")
        print("-s/--splitobs: create daily observation lists")
        print("-t/--test: test run, don't write data to files")
        print("datasetType: the dataset to be processed")

#------------------------------------------------------------------------------

def main(argv):
    dataType = ''
    splitObs = False
    testRun = False
    makeList = False
    removeDupl = False
    fileType = "fits"
    dataDir = defDataDir
    try:
        opts, args = getopt.getopt(argv[1:], "hd:f:lrst",
                                   ["help", "datadir:", "filetype:", "makelist",
                                    "splitobs", "remdupl", "test"])

    except getopt.GetoptError:
        # print help information and exit:
        print(argv)
        usage()
        raise SystemExit

    for o, a in opts:
        if o in ("-s", "--splitobs"):
            splitObs = True
        if o in ("-l", "--makelist"):
            makeList = True
        if o in ("-d", "--datadir"):
            dataDir = a
        if o in ("-f", "--filetype"):
            fileType = a        
        if o in ("-r", "--remdupl"):
            removeDupl = True
        if o in ("-h","--help"):
            usage()
            raise SystemExit
        if o in ("-t", "--test"):
            testRun = True

    if len(args) == 1:
        dataType = args[0]
    else:
        usage()
        raise SystemExit

    try:
        processData(dataType, makeList, dataDir, fileType, splitObs,
                    removeDupl, testRun)

    except Exception as e:
        print("An exception occurred:", str(e))
        with open(f'parse_to_ecsv_{dataType}_error.log', 'w') as f:
            with redirect_stdout(f):
                traceback.print_exc()

#------------------------------------------------------------------------------

if __name__ == "__main__":
    main(sys.argv)


#------------------------------------------------------------------------------
abcdefghilmqrstvwx

jknopuyz
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

##python src/autoIngest.py -s
python src/autoIngest.py -b
python src/autoIngest.py -w

python src/autoIngest.py -r c
python src/autoIngest.py -i c

python src/autoIngest.py -v

python src/autoIngest.py -r i
python src/autoIngest.py -i i

python src/autoIngest.py -c
python src/autoIngest.py -m

python src/autoIngest.py -r a
python src/autoIngest.py -i a





