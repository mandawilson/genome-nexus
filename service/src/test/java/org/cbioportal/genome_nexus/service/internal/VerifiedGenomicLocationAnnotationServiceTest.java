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
import static org.junit.Assert.fail;

import java.util.*;
import java.util.stream.Collectors;

import org.cbioportal.genome_nexus.model.VariantAnnotation;
import org.cbioportal.genome_nexus.service.GenomicLocationAnnotationService;
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
public class VerifiedGenomicLocationAnnotationServiceTest
{
    @InjectMocks
    private VerifiedGenomicLocationAnnotationServiceImpl verifiedGenomicLocationAnnotationServiceImpl;

    @Mock
    private GenomicLocationAnnotationService glVariantAnnotationService;

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
    public List<VariantTestCase> glSubstitutions = null;
    public List<VariantTestCase> glDeletions = null;
    public List<VariantTestCase> glInsertions = null;
    public List<VariantTestCase> glInsertionDeletions = null;
    // other types (duplications, no-change, methylation, specified alleles) are not handled by genome nexus and are not tested

    public String mockIsoformOverrideSource = "mockIsoformOverrideSource";
    public Map<String, String> mockTokenMap = new HashMap<String, String>();
    public List<String> mockFields = new ArrayList<String>();

    public void initializeTestCasesIfNeeded() {
        // VEP responses for these test cases were extraced from queries to http://genie.genomenexus.org/annotation/genomic/<variant>
        // these contain only the elements neccessary for testing the business logic in VerifiedGenomicLocationAnnotationService
        if (glSubstitutions == null) {
            glSubstitutions = new ArrayList<VariantTestCase>();
            glSubstitutions.add(new VariantTestCase("5,138163256,138163256,C,T", true, "C/T", true, "C/T", "valid substitution"));
            glSubstitutions.add(new VariantTestCase("5,138163256,138163256,A,T", false, null, false, null, "discrepant RefAllele"));
            glDeletions = new ArrayList<VariantTestCase>();
            glDeletions.add(new VariantTestCase("5,138163256,138163256,C,-", true, "C/-", true, "C/-", "1nt deletion with RefAllele"));
            glDeletions.add(new VariantTestCase("5,138163256,138163256,A,-", true, "C/-", false, null, "1nt deletion with discrepant RefAllele"));
            glDeletions.add(new VariantTestCase("5,138163255,138163256,TC,-", true, "TC/-", true, "TC/-", "2nt deletion with RefAllele"));
            glDeletions.add(new VariantTestCase("5,138163255,138163256,CC,-", true, "TC/-", false, null, "2nt deletion with discrepant RefAllele"));
            glDeletions.add(new VariantTestCase("5,138163255,138163256,CCCC,-", true, "TC/-", false, null, "2nt deletion with invalid RefAllele"));
            glInsertions = new ArrayList<VariantTestCase>();
            glInsertions.add(new VariantTestCase("5,138163255,138163256,-,T", true, "-/T", true, "-/T", "1nt insertion"));
            glInsertions.add(new VariantTestCase("5,138163255,138163256,-,TT", true, "-/TT", true, "-/TT", "2nt insertion"));
            glInsertions.add(new VariantTestCase("5,138163255,138163256,-,-", true, "-/-", true, "-/-", "insertion missing TumorSeqAllele")); // note : different result than ensembl
            glInsertionDeletions = new ArrayList<VariantTestCase>();
            glInsertionDeletions.add(new VariantTestCase("5,138163256,138163256,C,T", true, "C/T", true, "C/T", "1nt deletion with RefAllele, 1nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163256,138163256,A,T", true, "C/T", false, null, "1nt deletion with discrepant RefAllele, 1nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163256,138163256,C,TT", true, "C/TT", true, "C/TT", "1nt deletion with RefAllele, 2nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163256,138163256,A,TT", true, "C/TT", false, null, "1nt deletion with discrepant RefAllele, 2nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,TC,A", true, "TC/A", true, "TC/A", "2nt deletion with RefAllele, 1nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,TA,G", true, "TC/G", false, null, "2nt deletion with discrepant RefAllele, 1nt insertion"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,TC,TC", true, "TC/TC", true, null, "2nt deletion with RefAllele, 2nt insertion no change")); // note : different result than ensembl
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,TC,TT", true, "C/T", true, "C/T", "2nt deletion with RefAllele, 2nt insertion, partial change"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,TC,GG", true, "TC/GG", true, "TC/GG", "2nt deletion with RefAllele, 2nt insertion, full change"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,CC,TC", true, "TC/TC", false, null, "2nt deletion with discrepant RefAllele, 2nt insertion no change")); // note : different result than ensembl
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,CC,TT", true, "TC/TT", false, null, "2nt deletion with discrepant RefAllele, 2nt insertion, partial change"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,CC,GG", true, "TC/GG", false, null, "2nt deletion with discrepant RefAllele, 2nt insertion, full change"));
            glInsertionDeletions.add(new VariantTestCase("5,138163255,138163256,CCCC,TT", true, "TC/TT", false, null, "2nt deletion with invalid RefAllele, 2nt insertion"));
        }
    }

/*
    private void runTestSetByVariantString(List<VariantTestCase> variantTestCaseList, boolean supplyOtherParameters)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation testResponse = null;
            if (supplyOtherParameters) {
                testResponse = verifiedGenomicLocationAnnotationServiceImpl.getAnnotation(testCase.originalVariantQuery, mockIsoformOverrideSource, mockTokenMap, mockFields);
            } else {
                testResponse = verifiedGenomicLocationAnnotationServiceImpl.getAnnotation(testCase.originalVariantQuery);
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
*/
    private void runTestSetByVariantString(List<VariantTestCase> variantTestCaseList)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation testResponse = null;
            testResponse = verifiedGenomicLocationAnnotationServiceImpl.getAnnotation(testCase.originalVariantQuery);
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

/*
    private void runTestSetByVariantList(List<VariantTestCase> variantTestCaseList, boolean supplyOtherParameters)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        // make list of variant query strings and submit
        List<String> originalVariantQueries = variantTestCaseList.stream().map(t -> t.originalVariantQuery).collect(Collectors.toList());
        List<VariantAnnotation> variantResponse = null;
        if (supplyOtherParameters) {
            variantResponse = verifiedGenomicLocationAnnotationServiceImpl.getAnnotations(originalVariantQueries, mockIsoformOverrideSource, mockTokenMap, mockFields);
        } else {
            variantResponse = verifiedGenomicLocationAnnotationServiceImpl.getAnnotations(originalVariantQueries);
        }
        // check each element of response against expectations
        HashMap<String, VariantAnnotation> queryToResponse = new HashMap<String, VariantAnnotation>();
        for (VariantAnnotation responseElement : variantResponse) {
            if (queryToResponse.put(responseElement.getOriginalVariantQuery(), responseElement) != null) {
                fail("More than one response received for query string " + responseElement.getOriginalVariantQuery());
            }
        }
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation testResponse = queryToResponse.get(testCase.originalVariantQuery);
            if (testResponse == null) {
                fail("Response did not include record for query string " + testCase.originalVariantQuery);
            }
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
*/
    // Tests of the getAnnotation(String) function

    @Test
    public void getAnnotationForSubstitutionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glSubstitutions);
    }

    @Test
    public void getAnnotationForDeletionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glDeletions);
    }

    @Test
    public void getAnnotationForInsertionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glInsertions);
    }

    @Test
    public void getAnnotationForInsertionDeletionsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glInsertionDeletions);
    }
/*
    // Tests of the getAnnotation(String, overrideSource, tokenMap, fields) function

    @Test
    public void getAnnotationForSubstitutionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glSubstitutions, true);
    }

    @Test
    public void getAnnotationForDeletionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glDeletions, true);
    }

    @Test
    public void getAnnotationForInsertionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glInsertions, true);
    }

    @Test
    public void getAnnotationForInsertionDeletionsByVariantStringWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantString(glInsertionDeletions, true);
    }

    // Tests of the getAnnotations(List<String>) function

    @Test
    public void getAnnotationsForSubstitutionsByVariantList()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glSubstitutions, false);
    }

    @Test
    public void getAnnotationsForDeletionsByVariantList()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glDeletions, false);
    }

    @Test
    public void getAnnotationsForInsertionsByVariantList()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glInsertions, false);
    }

    @Test
    public void getAnnotationsForInsertionDeletionsByVariantList()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glInsertionDeletions, false);
    }

    // Tests of the getAnnotations(List<String>, overrideSource, tokenMap, fields) function

    @Test
    public void getAnnotationsForSubstitutionsByVariantListWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glSubstitutions, true);
    }

    @Test
    public void getAnnotationsForDeletionsByVariantListWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glDeletions, true);
    }

    @Test
    public void getAnnotationsForInsertionsByVariantListWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glInsertions, true);
    }

    @Test
    public void getAnnotationsForInsertionDeletionsByVariantListWithArgs()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        mockGenomicLocationAnnotationServiceMethods();
        runTestSetByVariantList(glInsertionDeletions, true);
    }
*/
    private void mockGenomicLocationAnnotationServiceMethodsForType(List<VariantTestCase> variantTestCaseList)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        ArrayList<String> queryList = new ArrayList<String>();
        ArrayList<VariantAnnotation> responseList = new ArrayList<VariantAnnotation>();
        for (VariantTestCase testCase : variantTestCaseList) {
            VariantAnnotation response = createAnnotation(
                    testCase.originalVariantQuery,
                    testCase.originalVariantQuery,
                    testCase.vepSuccessfullyAnnotated,
                    testCase.vepAlleleString);
            responseList.add(response);
            queryList.add(testCase.originalVariantQuery);
            Mockito.when(glVariantAnnotationService.getAnnotation(testCase.originalVariantQuery)).thenReturn(response);
            Mockito.when(glVariantAnnotationService.getAnnotation(testCase.originalVariantQuery, mockIsoformOverrideSource, mockTokenMap, mockFields)).thenReturn(response);
        }
        // order of response not guaranteed to match order of query - business logic should handle properly
        VariantAnnotation movedItem = responseList.remove(0);
        responseList.add(movedItem); // first element is now last
/*
        Mockito.when(glVariantAnnotationService.getAnnotations(queryList)).thenReturn(responseList);
        Mockito.when(glVariantAnnotationService.getAnnotations(queryList, mockIsoformOverrideSource, mockTokenMap, mockFields)).thenReturn(responseList);
*/
    }

    private void mockGenomicLocationAnnotationServiceMethods()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        initializeTestCasesIfNeeded();
        // mock methods in order to prevent hitting the live VEP web API
        mockGenomicLocationAnnotationServiceMethodsForType(glSubstitutions);
        mockGenomicLocationAnnotationServiceMethodsForType(glDeletions);
        mockGenomicLocationAnnotationServiceMethodsForType(glInsertions);
        mockGenomicLocationAnnotationServiceMethodsForType(glInsertionDeletions);
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
