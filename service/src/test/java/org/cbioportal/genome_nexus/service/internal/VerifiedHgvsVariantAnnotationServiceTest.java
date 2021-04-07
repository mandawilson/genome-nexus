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

    @Test
    public void getAnnotationByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        mockHgvsVariantAnnotationServiceMethods();
        String variant1 = "5:g.138163255_138163256delTCinsTT";
        VariantAnnotation test1 = variantAnnotationService.getAnnotation(variant1);
        assertEquals(test1.getOriginalVariantQuery(), variant1);
        assertTrue(test1.isSuccessfullyAnnotated());
        assertEquals(test1.getAlleleString(), "C/T");

        String variant2 = "5:g.138163255_138163256delCCinsTT";
        VariantAnnotation test2 = variantAnnotationService.getAnnotation(variant2);
        assertEquals(test2.getOriginalVariantQuery(), variant2);
        assertFalse(test2.isSuccessfullyAnnotated());

        String variant3 = "5:g.138163255_138163256delinsTT";
        VariantAnnotation test3 = variantAnnotationService.getAnnotation(variant3);
        assertEquals(test3.getOriginalVariantQuery(), variant3);
        assertTrue(test3.isSuccessfullyAnnotated());
        assertEquals(test3.getAlleleString(), "C/T");

    }

    @Test
    public void getAnnotationsByVariantString()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
    }


    private void mockHgvsVariantAnnotationServiceMethods()
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        // mock methods in order to prevent hitting the live VEP web API
        String variant1 = "5:g.138163255_138163256delTCinsTT";
        VariantAnnotation response1 = createAnnotation(variant1, variant1, true, "C/T");
        Mockito.when(hgvsVariantAnnotationService.getAnnotation(variant1)).thenReturn(response1);
        String variant2 = "5:g.138163255_138163256del";
        VariantAnnotation response2 = createAnnotation(variant2, variant2, true, "TC/-");
        Mockito.when(hgvsVariantAnnotationService.getAnnotation(variant2)).thenReturn(response2);
        String variant3 = "5:g.138163255_138163256delCCinsTT";
        VariantAnnotation response3 = createAnnotation(variant3, variant3, true, "C/T");
        Mockito.when(hgvsVariantAnnotationService.getAnnotation(variant3)).thenReturn(response3);
        String variant4 = "5:g.138163255_138163256delinsTT";
        VariantAnnotation response4 = createAnnotation(variant4, variant4, true, "C/T");
        Mockito.when(hgvsVariantAnnotationService.getAnnotation(variant4)).thenReturn(response4);
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
