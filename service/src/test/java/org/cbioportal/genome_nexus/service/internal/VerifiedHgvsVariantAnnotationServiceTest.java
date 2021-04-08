/*
 * Copyright (c) 2021 Memorial Sloan-Kettering Cancer Center.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS
 * FOR A PARTICULAR PURPOSE. The software and documentation provided hereunder
 * is on an "as is" basis, and Memorial Sloan-Kettering Cancer Center has no
 * obligations to provide maintenance, support, updates, enhancements or
 * modifications. In no event shall Memorial Sloan-Kettering Cancer Center be
 * liable to any party for direct, indirect, special, incidental or
 * consequential damages, including lost profits, arising out of the use of this
 * software and its documentation, even if Memorial Sloan-Kettering Cancer
 * Center has been advised of the possibility of such damage.
 */

/*
 * This file is part of cBioPortal.
 *
 * cBioPortal is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package org.cbioportal.genome_nexus.service.internal;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import java.util.*;

import org.cbioportal.genome_nexus.model.VariantAnnotation;
import org.cbioportal.genome_nexus.service.exception.VariantAnnotationNotFoundException;
import org.cbioportal.genome_nexus.service.exception.VariantAnnotationWebServiceException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.Spy;
import org.mockito.junit.MockitoJUnitRunner;

@RunWith(MockitoJUnitRunner.Silent.class)
public class VerifiedHgvsVariantAnnotationServiceTest
{
    @InjectMocks
    private VerifiedHgvsVariantAnnotationService variantAnnotationService;

    @Mock
    private HgvsVariantAnnotationService hgvsVariantAnnotationService;

    class VariantTestCase {
        public String originalVariantQuery;
        public boolean vepSuccessfullyAnnotated;
        public String vepAlleleString;
        public boolean expectedGnSuccessfullyAnnotated;
        public String expectedGnAlleleString;
        public String description;
        public VariantTestCase(
                String originalVariantQuery,
                boolean vepSuccessfullyAnnotated,
                String vepAlleleString,
                boolean expectedGnSuccessfullyAnnotated,
                String expectedGnAlleleString,
                String description) {
            this.originalVariantQuery =  originalVariantQuery;
            this.vepSuccessfullyAnnotated = vepSuccessfullyAnnotated;
            this.vepAlleleString = vepAlleleString;
            this.expectedGnSuccessfullyAnnotated = expectedGnSuccessfullyAnnotated;
            this.expectedGnAlleleString = expectedGnAlleleString;
            this.description = description;
        }
    }

    // test cases
    public List<VariantTestCase> hgvsSubstitutions = null;
    public List<VariantTestCase> hgvsDeletions = null;
    public List<VariantTestCase> hgvsInsertions = null;
    public List<VariantTestCase> hgvsInsertionDeletions = null;
    public List<VariantTestCase> hgvsInversions = null; // test results of an unsupported format query
    // other types (duplications, no-change, methylation, specified alleles) are not handled by genome nexus and are not tested

    public String mockIsoformOverrideSource = "mockIsoformOverrideSource";
    public Map<String, String> mockTokenMap = new HashMap<String, String>();
    public List<String> mockFields = new ArrayList<String>();

    public void initializeTestCasesIfNeeded() {
        // VEP responses for these test cases are extraced from queries to http://grch37.rest.ensembl.org/vep/human/hgvs/<variant>
        if (hgvsSubstitutions == null) {
            hgvsSubstitutions = new ArrayList<VariantTestCase>();
            hgvsSubstitutions.add(new VariantTestCase("5:g.138163256C>T", true, "C/T", true, "C/T", "valid substitution"));
            hgvsSubstitutions.add(new VariantTestCase("5:g.138163256A>T", false, null, false, null, "discrepant RefAllele"));
            hgvsSubstitutions.add(new VariantTestCase("5:g.138163256>T", false, null, false, null, "missing RefAllele"));
            hgvsDeletions = new ArrayList<VariantTestCase>();
            hgvsDeletions.add(new VariantTestCase("5:g.138163256delC", true, "C/-", true, "C/-", "1nt deletion with RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163256delA", true, "C/-", false, null, "1nt deletion with discrepant RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163256del", true, "C/-", true, "C/-", "1nt deletion missing RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163255_138163256delTC", true, "TC/-", true, "TC/-", "2nt deletion with RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163255_138163256delCC", true, "TC/-", false, null, "2nt deletion with discrepant RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163255_138163256delCCCC", true, "TC/-", false, null, "2nt deletion with invalid RefAllele"));
            hgvsDeletions.add(new VariantTestCase("5:g.138163255_138163256del", true, "TC/-", true, "TC/-", "2nt deletion missing RefAllele"));
            hgvsInsertions = new ArrayList<VariantTestCase>();
            hgvsInsertions.add(new VariantTestCase("5:g.138163255_138163256insT", true, "-/T", true, "-/T", "1nt insertion"));
            hgvsInsertions.add(new VariantTestCase("5:g.138163255_138163256insTT", true, "-/TT", true, "-/TT", "2nt insertion"));
            hgvsInsertions.add(new VariantTestCase("5:g.138163255_138163256ins", false, null, false, null, "insertion missing TumorSeqAllele"));
            hgvsInsertionDeletions = new ArrayList<VariantTestCase>();
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delinsT", true, "C/T", true, "C/T", "1nt deletion without RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delCinsT", true, "C/T", true, "C/T", "1nt deletion with RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delAinsT", true, "C/T", false, null, "1nt deletion with discrepant RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delinsC", false, null, false, null, "1nt deletion without RefAllele, 1nt insertion no change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delinsTT", true, "C/TT", true, "C/TT", "1nt deletion without RefAllele, 2nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delCinsTT", true, "C/TT", true, "C/TT", "1nt deletion with RefAllele, 2nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163256delAinsTT", true, "C/TT", false, null, "1nt deletion with discrepant RefAllele, 2nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delinsG", true, "TC/G", true, "TC/G", "2nt deletion without RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delTCinsA", true, "TC/A", true, "TC/A", "2nt deletion with RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delTAinsG", true, "TC/G", false, null, "2nt deletion with discrepant RefAllele, 1nt insertion"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delTCinsTC", false, null, false, null, "2nt deletion with RefAllele, 2nt insertion no change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delTCinsTT", true, "C/T", true, "C/T", "2nt deletion with RefAllele, 2nt insertion, partial change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delTCinsGG", true, "TC/GG", true, "TC/GG", "2nt deletion with RefAllele, 2nt insertion, full change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delCCinsTC", false, null, false, null, "2nt deletion with discrepant RefAllele, 2nt insertion no change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delCCinsTT", true, "C/T", false, null, "2nt deletion with discrepant RefAllele, 2nt insertion, partial change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delCCinsGG", true, "TC/GG", false, null, "2nt deletion with discrepant RefAllele, 2nt insertion, full change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delinsTC", false, null, false, null, "2nt deletion without RefAllele, 2nt insertion no change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delinsTT", true, "C/T", true, "C/T", "2nt deletion without RefAllele, 2nt insertion, partial change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delinsGG", true, "TC/GG", true, "TC/GG", "2nt deletion without RefAllele, 2nt insertion, full change"));
            hgvsInsertionDeletions.add(new VariantTestCase("5:g.138163255_138163256delCCCCinsTT", true, "C/T", false, null, "2nt deletion with invalid RefAllele, 2nt insertion"));
            hgvsInversions = new ArrayList<VariantTestCase>();
            hgvsInversions.add(new VariantTestCase("5:g.138163255_138163256inv", true, "TC/GA", true, "TC/GA", "inversions not supported - but will run as passthrough"));
            hgvsInversions.add(new VariantTestCase("5:g.138163255_138163256invTC", false, null, false, null, "inversion format does not allow specification of RefAllele"));
        }
    }

    private void runTestSetByVariantString(List<VariantTestCase> variantTestCaseList, boolean supplyOtherParameters)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation testResponse = null;
            if (supplyOtherParameters) {
                testResponse = variantAnnotationService.getAnnotation(testCase.originalVariantQuery, mockIsoformOverrideSource, mockTokenMap, mockFields);
            } else {
                testResponse = variantAnnotationService.getAnnotation(testCase.originalVariantQuery);
            }
            assertEquals(testCase.originalVariantQuery + " : response query field does not match request query string", testCase.originalVariantQuery, testResponse.getOriginalVariantQuery());
            if (testCase.expectedGnSuccessfullyAnnotated) {
                assertTrue(testCase.originalVariantQuery + " : expected successful annotation", testResponse.isSuccessfullyAnnotated());
            } else {
                assertFalse(testCase.originalVariantQuery + " : expected failed annotation", testResponse.isSuccessfullyAnnotated());
            }
            if (testResponse.isSuccessfullyAnnotated()) {
                assertEquals(testCase.originalVariantQuery + " : Variant Allele comparison", testCase.expectedGnAlleleString, testResponse.getAlleleString());
            }
        }
    }

    // Tests of the getVarant(String) function

    @Test
    public void getAnnotationForSubstitutionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsSubstitutions, false);
    }

    @Test
    public void getAnnotationForDeletionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsDeletions, false);
    }

    @Test
    public void getAnnotationForInsertionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInsertions, false);
    }

    @Test
    public void getAnnotationForInsertionDeletionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInsertionDeletions, false);
    }

    @Test
    public void getAnnotationForInversionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInversions, false);
    }

    // Tests of the getVarant(String, overrideSource, tokenMap, fields) function

    @Test
    public void getAnnotationForSubstitutionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsSubstitutions, true);
    }

    @Test
    public void getAnnotationForDeletionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsDeletions, true);
    }

    @Test
    public void getAnnotationForInsertionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInsertions, true);
    }

    @Test
    public void getAnnotationForInsertionDeletionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInsertionDeletions, true);
    }

    @Test
    public void getAnnotationForInversionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockHgvsVariantAnnotationServiceMethods();
        runTestSetByVariantString(hgvsInversions, true);
    }

    private void mockHgvsVariantAnnotationServiceMethodsForType(List<VariantTestCase> variantTestCaseList)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation response = createAnnotation(
                    testCase.originalVariantQuery,
                    testCase.originalVariantQuery,
                    testCase.vepSuccessfullyAnnotated,
                    testCase.vepAlleleString);
            Mockito.when(hgvsVariantAnnotationService.getAnnotation(testCase.originalVariantQuery)).thenReturn(response);
            Mockito.when(hgvsVariantAnnotationService.getAnnotation(testCase.originalVariantQuery, mockIsoformOverrideSource, mockTokenMap, mockFields)).thenReturn(response);
        }
    }

    private void mockHgvsVariantAnnotationServiceMethods()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        // mock methods in order to prevent hitting the live VEP web API
        mockHgvsVariantAnnotationServiceMethodsForType(hgvsSubstitutions);
        mockHgvsVariantAnnotationServiceMethodsForType(hgvsDeletions);
        mockHgvsVariantAnnotationServiceMethodsForType(hgvsInsertions);
        mockHgvsVariantAnnotationServiceMethodsForType(hgvsInsertionDeletions);
        mockHgvsVariantAnnotationServiceMethodsForType(hgvsInversions);
    }

    private VariantAnnotation createAnnotation(String originalVariantQuery, String variant, boolean successfullyAnnotated, String alleleString)
    {
        VariantAnnotation annotation = new VariantAnnotation();
        annotation.setOriginalVariantQuery(originalVariantQuery);
        annotation.setVariant(variant);
        annotation.setAlleleString(alleleString);
        annotation.setSuccessfullyAnnotated(successfullyAnnotated);
        return annotation;
    }

}
