from unittest.mock import MagicMock, patch

import pytest

from src.defs.auto_acmg import (
    AutoACMGCriteria,
    AutoACMGPrediction,
    AutoACMGStrength,
    AutoACMGStrucVarData,
    GenomicStrand,
)
from src.defs.auto_pvs1 import PVS1Prediction, PVS1PredictionStrucVarPath
from src.defs.exceptions import AlgorithmError, MissingDataError
from src.defs.genome_builds import GenomeRelease
from src.defs.mehari import Exon
from src.defs.strucvar import StrucVar, StrucVarType
from src.strucvar.auto_pvs1 import AutoPVS1, StrucVarHelper


@pytest.fixture
def strucvar_helper():
    return StrucVarHelper()


@pytest.fixture
def strucvar():
    return StrucVar(
        sv_type=StrucVarType.DEL,
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=100,
        stop=200,
    )


@pytest.fixture
def exons():
    return [MagicMock(altStartI=50, altEndI=150), MagicMock(altStartI=160, altEndI=210)]


# ========== StrucVarHelper ============

# --------- full_gene_del ---------


def test_full_gene_del_true(strucvar_helper, strucvar, exons):
    """Test if full_gene_del correctly identifies a full gene deletion."""
    # Modify strucvar to fully encompass the gene defined by exons
    strucvar.start = 45
    strucvar.stop = 215
    assert (
        strucvar_helper.full_gene_del(strucvar, exons) is True
    ), "Should recognize a full gene deletion"


def test_full_gene_del_false(strucvar_helper, strucvar, exons):
    """Test if full_gene_del correctly identifies when it's not a full gene deletion."""
    # Modify strucvar to not fully encompass the gene defined by exons
    strucvar.start = 60
    strucvar.stop = 205
    assert (
        strucvar_helper.full_gene_del(strucvar, exons) is False
    ), "Should recognize not a full gene deletion"


def test_full_gene_del_missing_exons(strucvar_helper, strucvar):
    """Test full_gene_del with no exons provided."""
    with pytest.raises(MissingDataError):
        strucvar_helper.full_gene_del(strucvar, [])


def test_full_gene_del_exon_boundaries(strucvar_helper, strucvar):
    """Test the edge cases of exon boundaries."""
    exons = [MagicMock(altStartI=100, altEndI=200)]
    strucvar.start = 100
    strucvar.stop = 200
    assert (
        strucvar_helper.full_gene_del(strucvar, exons) is True
    ), "Should handle edge case boundaries correctly"

    strucvar.start = 99
    strucvar.stop = 201
    assert (
        strucvar_helper.full_gene_del(strucvar, exons) is True
    ), "Should handle gene deletion touching boundaries"

    strucvar.start = 101
    strucvar.stop = 199
    assert (
        strucvar_helper.full_gene_del(strucvar, exons) is False
    ), "Should not consider partial deletions"


# --------- del_disrupt_rf ---------


def test_del_disrupt_rf_full_exon_deletion(strucvar_helper, strucvar):
    """Test if del_disrupt_rf correctly identifies a full exon deletion."""
    exons = [
        MagicMock(altStartI=50, altEndI=100),
        MagicMock(altStartI=150, altEndI=200),
        MagicMock(altStartI=250, altEndI=300),
    ]
    strucvar.start = 140
    strucvar.stop = 210
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is False
    ), "Full exon deletion should not disrupt reading frame"


def test_del_disrupt_rf_partial_exon_deletion_plus_strand(strucvar_helper, strucvar):
    """Test if del_disrupt_rf correctly identifies a partial exon deletion on plus strand."""
    exons = [
        MagicMock(altStartI=50, altEndI=80),
        MagicMock(altStartI=100, altEndI=184),
        MagicMock(altStartI=200, altEndI=300),
        MagicMock(altStartI=320, altEndI=350),
    ]

    # Deletion of 3 bases at the start of the first exon and the entire second exon
    strucvar.start = 52
    strucvar.stop = 197
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is False
    ), "Deletion of 3 bases at exon end and a full exon should not disrupt reading frame"

    # Deletion of 4 bases at the start of the first exon and the entire second exon
    strucvar.start = 53
    strucvar.stop = 197
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is True
    ), "Deletion of 4 bases at exon end and a full exon should disrupt reading frame"

    # Deletion of 3 bases at the start of the second exon and the entire first exon
    strucvar.start = 30
    strucvar.stop = 102
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is False
    ), "Deletion of 3 bases at exon start and a full exon should not disrupt reading frame"

    # Deletion of 4 bases at the start of the second exon and the entire first exon
    strucvar.start = 30
    strucvar.stop = 103
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is True
    ), "Deletion of 4 bases at exon start and a full exon should disrupt reading frame"


def test_del_disrupt_rf_partial_exon_deletion_minus_strand(strucvar_helper, strucvar):
    """Test if del_disrupt_rf correctly identifies a partial exon deletion on minus strand."""
    exons = [MagicMock(altStartI=100, altEndI=200), MagicMock(altStartI=250, altEndI=300)]

    # Deletion of 3 bases at the start of the exon
    strucvar.start = 100
    strucvar.stop = 298
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Minus) is False
    ), "Deletion of 3 bases at exon start should not disrupt reading frame on minus strand"

    # Deletion of 4 bases at the start of the exon
    strucvar.start = 100
    strucvar.stop = 297
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Minus) is True
    ), "Deletion of 4 bases at exon start should disrupt reading frame on minus strand"

    # Deletion of 3 bases at the end of the exon
    strucvar.start = 198
    strucvar.stop = 300
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Minus) is False
    ), "Deletion of 3 bases at exon end should not disrupt reading frame on minus strand"

    # Deletion of 4 bases at the end of the exon
    strucvar.start = 197
    strucvar.stop = 300
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Minus) is True
    ), "Deletion of 4 bases at exon end should disrupt reading frame on minus strand"


def test_del_disrupt_rf_multiple_exons(strucvar_helper, strucvar):
    """Test if del_disrupt_rf correctly handles deletions spanning multiple exons."""
    exons = [
        MagicMock(altStartI=100, altEndI=200),
        MagicMock(altStartI=300, altEndI=400),
        MagicMock(altStartI=500, altEndI=600),
    ]
    strucvar.start = 150
    strucvar.stop = 350
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is False
    ), "Deletion of full exons should not disrupt reading frame"


def test_del_disrupt_rf_intronic_deletion(strucvar_helper, strucvar):
    """Test if del_disrupt_rf correctly handles intronic deletions."""
    exons = [
        MagicMock(altStartI=100, altEndI=200),
        MagicMock(altStartI=300, altEndI=400),
    ]
    strucvar.start = 1
    strucvar.stop = 299
    assert (
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus) is False
    ), "Intronic deletion should not disrupt reading frame"


def test_del_disrupt_rf_no_affected_exons(strucvar_helper, strucvar):
    """Test if del_disrupt_rf raises an error when no exons are affected."""
    exons = [MagicMock(altStartI=100, altEndI=200)]
    strucvar.start = 50
    strucvar.stop = 99
    with pytest.raises(AlgorithmError):
        strucvar_helper.del_disrupt_rf(strucvar, exons, GenomicStrand.Plus)


# ========== AutoPVS1 ============


@pytest.fixture
def auto_pvs1():
    return AutoPVS1()


@pytest.fixture
def strucvar_del():
    return StrucVar(
        sv_type=StrucVarType.DEL,
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=100,
        stop=200,
    )


@pytest.fixture
def strucvar_dup():
    return StrucVar(
        sv_type=StrucVarType.DUP,
        genome_release=GenomeRelease.GRCh38,
        chrom="1",
        start=100,
        stop=200,
    )


@pytest.fixture
def var_data():
    data = AutoACMGStrucVarData()
    data.exons = [MagicMock(altStartI=50, altEndI=150), MagicMock(altStartI=160, altEndI=210)]
    return data


@patch.object(AutoPVS1, "full_gene_del", return_value=True)
def test_verify_pvs1_full_gene_deletion(mock_full_gene_del, auto_pvs1, strucvar_del, var_data):
    """Test PVS1 when a full gene deletion is detected."""
    prediction, path, comment = auto_pvs1.verify_pvs1(strucvar_del, var_data)
    assert prediction == PVS1Prediction.PVS1, "Should predict PVS1 for full gene deletions"


@patch.object(AutoPVS1, "del_disrupt_rf", return_value=True)
@patch.object(AutoPVS1, "undergo_nmd", return_value=False)
@patch.object(AutoPVS1, "crit4prot_func", return_value=True)
def test_verify_pvs1_deletion_disrupts_rf(
    mock_crit4prot_func, mock_undergo_nmd, mock_del_disrupt_rf, auto_pvs1, strucvar_del, var_data
):
    """Test PVS1 when a deletion disrupts the reading frame but does not undergo NMD."""
    prediction, path, comment = auto_pvs1.verify_pvs1(strucvar_del, var_data)
    assert (
        prediction == PVS1Prediction.PVS1_Strong
    ), "Should predict PVS1 Strong for disruptive deletions not undergoing NMD"


@patch.object(AutoPVS1, "proven_in_tandem", return_value=True)
@patch.object(AutoPVS1, "dup_disrupt_rf", return_value=True)
@patch.object(AutoPVS1, "undergo_nmd", return_value=True)
def test_verify_pvs1_duplication_in_tandem(
    mock_undergo_nmd, mock_dup_disrupt_rf, mock_proven_in_tandem, auto_pvs1, strucvar_dup, var_data
):
    """Test PVS1 when a duplication is in tandem and disrupts the reading frame."""
    prediction, path, comment = auto_pvs1.verify_pvs1(strucvar_dup, var_data)
    assert (
        prediction == PVS1Prediction.PVS1
    ), "Should predict PVS1 for tandem duplications disrupting RF and undergoing NMD"


@patch.object(
    AutoPVS1,
    "verify_pvs1",
    return_value=(
        PVS1Prediction.PVS1,
        PVS1PredictionStrucVarPath.DEL1,
        "Full gene deletion identified",
    ),
)
def test_predict_pvs1(mock_verify_pvs1, auto_pvs1, strucvar_del, var_data):
    """Test the predict_pvs1 method for structural variants."""
    criteria = auto_pvs1.predict_pvs1(strucvar_del, var_data)
    assert isinstance(criteria, AutoACMGCriteria), "Should return an instance of AutoACMGCriteria"
    assert criteria.prediction == AutoACMGPrediction.Met, "Prediction should be Met"
    assert (
        criteria.strength == AutoACMGStrength.PathogenicVeryStrong
    ), "Strength should be Very Strong"
    assert (
        "Full gene deletion identified" in criteria.summary
    ), "Summary should include comment from verification"
