/*
 * Copyright (c) 2021
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 */

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimMetaBall.h>
#include <GU/GU_PrimPart.h>
#include <GA/GA_Primitive.h>
#include <UT/UT_VectorTypes.h>
#include <GU/GU_PrimMetaBall.h>
#include <GU/GU_PrimSphere.h>
#include <stddef.h>
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int getline(char **lineptr, size_t *n, FILE *stream)
{
static char line[256];
char *ptr;
unsigned int len;

   if (lineptr == NULL || n == NULL)
   {
      errno = EINVAL;
      return -1;
   }

   if (ferror (stream)) {
       printf("failure 2\n");
      return -1;
   }

   if (feof(stream)) {
       printf("failure 3\n");
      return -1;
   }
     
   fgets(line,256,stream);

   ptr = strchr(line,'\n');   
   if (ptr)
      *ptr = '\0';

   len = strlen(line);
   
   if ((len+1) < 256)
   {
      ptr = (char*)realloc(*lineptr, 256);
      if (ptr == NULL) {
          printf("failure 4\n");
         return(-1);
      }
      *lineptr = ptr;
      *n = 256;
   }

   strcpy(*lineptr,line); 
   return(len);
}

namespace HDK_Sample {

static float
densityFunction(const UT_Vector3 &P, void *data)
{
    // Return the signed distance to the unit sphere
    return 1 - P.length();
}

}

void
sopCreateBall(GU_Detail *to, float x, float y, float z)
{
    if(true) {
        GA_Offset ptoff = to->appendPointOffset();
        to->setPos3(ptoff, x, y, z);
        return;
    }
    if(true) {
        GU_PrimMetaBallParms ball(to);
        ball.weight = 200;
        // Position and scale the collision sphere to account for the
        // deformation animation.
        ball.xform.scale(0.005, 0.005, 0.005);
        ball.xform.translate(x, y, z);
        GU_PrimMetaBall::build(ball, "blinn");
    } else {
        GU_PrimSphereParms ball(to);
        //ball.weight = 200;
        // Position and scale the collision sphere to account for the
        // deformation animation.
        ball.xform.scale(0.005, 0.005, 0.005);
        ball.xform.translate(x, y, z);
        GU_PrimSphere::build(ball);
    }
}

int
main(int argc, char *argv[])
{
    assert(argc == 3);
    char *fname = argv[1];
    FILE * fp;
    fp = fopen(fname, "r");
    if(fp == NULL) {
        printf("Failed to open %s\n", fname);
        exit(1);
    }
    char *out_fname = (char*)malloc(strlen(argv[2]) + 20);
    int i = 0;
    char *line = NULL; 
    size_t len = 0;
    int read;
    GU_Detail gdp;

    getline(&line, &len, fp);
    int nparticles;
    sscanf(line, "%d\n", &nparticles);
    float *particles = (float*)malloc(sizeof(float)*nparticles*3);
    GEO_PrimParticle *sys;
    sys = GU_PrimParticle::build(&gdp, nparticles);
    UT_Matrix3F scale;
    scale.identity();
    scale.scale(0.05, 0.05, 0.05);
    /*GEO_PrimSphere **balls = (GEO_PrimSphere**)malloc(sizeof(void *) * nparticles);
    for(int i = 0; i < nparticles; i++) {//nparticles; i++) {
        balls[i] = (GEO_PrimSphere *)gdp.appendPrimitive(GA_PRIMSPHERE);
        //balls[i]->setWeight(1);
        //balls[i]->setTransform(scale);
        //ball->setPos3(UT_Vector3())
    }*/
    sys->getRenderAttribs().setType(GEO_PARTICLE_SPHERE);
    sys->getRenderAttribs().setMotionBlur(1);
    sys->getRenderAttribs().setBlurTime(0.03);
    sys->getRenderAttribs().setSize(1);
    //GA_Primitive::const_iterator it;
    //sys->beginVertex(it);
    int pidx = 0;
    int every = 1;
    while((read = getline(&line, &len, fp)) != -1) {
        //printf("Read line %s\n", line);
        if(strncmp("---DONE---", line, 10) == 0) {
            
            if(i % every == 0) {
                assert(pidx == nparticles);
                sprintf(out_fname, "%s\\frame%d.bgeo", argv[2], i);
            
                printf("Writing out to %s\n", out_fname);
            
                gdp.save(out_fname, NULL);
                gdp.clearAndDestroy();
            }
            i++;
            pidx = 0;
            //sys->beginVertex(it);
            continue;
        }
        if(i % every != 0)
            continue;
        float x, y, z;
        sscanf(line, "%f,%f,%f\n", &x, &y, &z);
        particles[pidx*3] = x;
        particles[pidx*3+1] = y;
        particles[pidx*3+2] = z;
        //printf("Read line %s\n", line);
        //gdp.setPos3(sys->getPointOffset(pidx), x, y, z);
        UT_Matrix4F scale;
        scale.identity();
        scale.scale(0.05, 0.05, 0.05);
        scale.translate(x, y, z);
        //balls[pidx]->setTransform4(scale);
        //balls[pidx]->setPos3(UT_Vector3(x, y, z));
        sopCreateBall(&gdp, x, y, z);
        //printf("Read line %s\n", line);
        //it.advance();
        pidx++;
    }
    return 0;
    /*GU_Detail		gdp;
    UT_BoundingBox	bounds;

    // Evaluate the iso-surface inside this bounding box
    bounds.setBounds(-1, -1, -1, 1, 1, 1);

    // Create an iso-surface
    gdp.polyIsoSurface(HDK_Sample::densityFunction, NULL, bounds, 20, 20, 20);

    // Save to sphere.bgeo
    gdp.save("sphere.bgeo", NULL);*/

    return 0;
}
