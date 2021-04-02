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

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.cbioportal.genome_nexus.service.*;
import org.cbioportal.genome_nexus.service.internal.HgvsVariantAnnotationService;
import org.cbioportal.genome_nexus.service.exception.VariantAnnotationNotFoundException;
import org.cbioportal.genome_nexus.service.exception.VariantAnnotationWebServiceException;
import org.cbioportal.genome_nexus.model.VariantAnnotation;
import org.cbioportal.genome_nexus.component.annotation.NotationConverter;
import org.cbioportal.genome_nexus.service.cached.CachedVariantAnnotationFetcher;
import org.springframework.context.annotation.Lazy;
import org.springframework.stereotype.Service;
import org.springframework.beans.factory.annotation.*;
import java.util.*;

@Service
public class VerifiedHgvsVariantAnnotationService implements VariantAnnotationService
{
    private static final Log LOG = LogFactory.getLog(VerifiedHgvsVariantAnnotationService.class);
    private final HgvsVariantAnnotationService hgvsVariantAnnotationService;

    @Autowired
    public VerifiedHgvsVariantAnnotationService(
        HgvsVariantAnnotationService hgvsVariantAnnotationService)
    {
        this.hgvsVariantAnnotationService = hgvsVariantAnnotationService;
    }

    @Override
    public VariantAnnotation getAnnotation(String variant)
        throws VariantAnnotationNotFoundException, VariantAnnotationWebServiceException
    {
        return hgvsVariantAnnotationService.getAnnotation(variant);
    }

    @Override
    public List<VariantAnnotation> getAnnotations(List<String> variants)
    {
        return hgvsVariantAnnotationService.getAnnotations(variants);
    }

    @Override
    public VariantAnnotation getAnnotation(String variant, String isoformOverrideSource, Map<String, String> token, List<String> fields)
        throws VariantAnnotationWebServiceException, VariantAnnotationNotFoundException
    {
        return hgvsVariantAnnotationService.getAnnotation(variant, isoformOverrideSource, token, fields);
    }

    @Override
    public List<VariantAnnotation> getAnnotations(List<String> variants, String isoformOverrideSource, Map<String, String> token, List<String> fields)
    {
        return hgvsVariantAnnotationService.getAnnotations(variants, isoformOverrideSource, token, fields);
    }
}
