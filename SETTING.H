/****************************************************************************/
/*
 * Function name: newProject
 * Description: creat a project file
 * Parameters:
 *			 proj_path
 * Return:
 *			 1				Successful
 *			 0				Fail
 *			 -1             No project created
 */
/* If compute radiation power, define _POYNTING to compile the module */
#define _POYNTING
/* If compute SAR, define _SAR to compile the module */
#define _SAR
int newProject(void)
{
	int i;
	FILE* fp_proj;

/* Begin setting*/
/****************************************************************************/
/*  1. Project                                                              */
/*  If creat a new project, set isNewProj = 1                               */
/*  Else(Use an existing project) set isNewProj = 0, completed block 1, then*/
/*  compile.(skip others block)                                             */
/****************************************************************************/

    isNewProj = 0;

    char set_path_proj[] = "F:\\work\\TestData\\";
    char set_proj_name[] = "MedFDTD_Proj";
    char set_path_log[] = "F:\\work\\TestData\\";
    char set_path_save[] = "F:\\work\\TestData\\";

/***END OF BLOCK 1***********************************************************/

/****************************************************************************/
/*  2. Mesh                                                                 */
/****************************************************************************/

	/* Size of computational grids(x, y, z) */
	_spaceX=20;
	_spaceY=20;
	_spaceZ=90;
    /* Space from absorbing boundary to computational grid */
    padding=0;
    /* Spatial resolution(Meter) */
    dx = 0.001;
    dy = 0.001;
    dz = 0.001;

/***END OF BLOCK 2***********************************************************/

/****************************************************************************/
/*  3. Simulation time                                                      */
/****************************************************************************/

    /* Max time step */
    nMax=2000;
    /* Time step */
    dt = 1.0 / (C * sqrt(1.0/(dx * dx) + 1.0/(dy * dy) + 1.0/(dz * dz)));/* The maximum value of dt */

/***END OF BLOCK 3***********************************************************/

/****************************************************************************/
/*  4. Absorbing boundary                                                   */
/****************************************************************************/

    /* Type of absorbing boundary : 0-->PEC
                                    1-->PML
                                    2-->Mur2(only available in serial)
    */
	abcNo = 1;
	thicknessOfPml=7;

/***END OF BLOCK 4***********************************************************/

/****************************************************************************/
/*  5. Power source                                                         */
/****************************************************************************/

    /************************************************************************/
    /*  5.1 Source type                                                     */
    /************************************************************************/

    /* Set source type: 0-->Point source
                        1-->Plane source(if choose 1, please go to block 5.3)
                        2-->Dipole source
    */
    sourceType=0;
	
    /***END OF BLOCK 5.1*****************************************************/

    /************************************************************************/
    /*  5.2 Port                                                            */
    /************************************************************************/

    /* Location of source: range 1 to _spaceX(_spaceY or _spaceZ) in Block 2*/
	_isource=11;
	_jsource=11;
	_ksource=45;
    /* Voltage direction, 'x'-->Ex, 'y'-->Ey, 'z'-->Ez */
	port='z';

    /***END OF BLOCK 5.2*****************************************************/

    /************************************************************************/
    /*  5.3 Wave form                                                       */
    /************************************************************************/

    /* Wave form: -1-->User defined : pathSRC_set (Available parameters)
                   0-->Sine : freq (Available parameters)
                   1-->Gauss : t0, pulse_width
                   2-->raised cosine : pulse_width;
	 		       3-->differential Gauss : t0, pulse_width
                   4-->3-truncated cosine : pulse_width; 
	 		       5-->modulation Gauss : freq, t0, pulse_width.
	*/
	waveForm=0;
	amp=1;

	freq=1.8e9;
	t0 = 150;
	pulse_width = 120;

    char set_pathSRC[] = "C:\\source.txt";

    /***END OF BLOCK 5.3*****************************************************/

/***END OF BLOCK 5***********************************************************/

/****************************************************************************/
/*  6. Model                                                                */
/****************************************************************************/

    /************************************************************************/
    /*  6.1 Import model                                                    */
    /************************************************************************/
    /* Path to import model */
    char set_path_data[] = "F:\\ModelData";
    char set_model_name[] = "model.txt";/* Order: z(outer)<--y<--x(inner) */
    char set_media_name[] = "media.txt";/* Form: id sigma epsilon rho\n */
    /* Number of media type in media file */
	mediaNum = 0;
    /***END OF BLOCK 6.1*****************************************************/

    /************************************************************************/
    /*  6.2 Build object                                                    */
    /*  Cover the import model if overlap in space                          */
    /************************************************************************/
    object_num = 1;
    int object_data[][8] = 
    {
        13, 2, 5, 4, 6, 90, 0
    };
    double new_media[][3] = 
    {
        0.1, 81, 1e6
    };
    /***END OF BLOCK 6.2*****************************************************/

    /************************************************************************/
    /*  6.3 Antenna                                                         */
    /************************************************************************/
    int antenna_amount = 1;
	int antenna_direction = 3;
	int antenna_feed_x = _isource;
	int antenna_feed_y = _jsource;
	int antenna_feed_z = _ksource;
	int antenna_impedance = 0;
	int antenna_length_high = 41;
	int antenna_length_low  = 41;
    /***END OF BLOCK 6.3*****************************************************/

/***END OF BLOCK 6***********************************************************/

/****************************************************************************/
/*  7. Save Field Snapshots Setting:                                        */
/****************************************************************************/

	save_plane_amount = 3;
    int save_plane_parameters[][5] = 
    {/* Input saving parameters
        Form: Begining Time Step, Ending Time Step, Increment, Plane, Slice.
        Plane No.: 1-->xy
                   2-->yz
                   3-->xz
     */
        1, nMax, 1, 1, _ksource,
        1, nMax, 1, 2, _isource,
        1, nMax, 1, 3, _jsource
    };

/***END OF BLOCK 7***********************************************************/

#ifdef _SAR
/****************************************************************************/
/*  8. SAR Setting:                                                         */
/****************************************************************************/

    /************************************************************************/
    /*  8.1 Local SAR Setting:                                              */
    /*  If need compute mass averaged sar, please go to block 7.2           */
    /************************************************************************/

    save_localSAR_amount = 3;
    int save_sar_parameters[][3] = 
    {/* Input saving parameters
        Form: Begining Time Step, Ending Time Step, Slice.
     */
        nMax-(int)(1/freq/dt)+1, nMax, _ksource-1,/* (int)(1/freq/dt): a period */
        nMax-(int)(1/freq/dt)+1, nMax, _ksource,
        nMax-(int)(1/freq/dt)+1, nMax, _ksource+1
    };

    /***END OF BLOCK 8.1*****************************************************/

    /************************************************************************/
    /*  8.2 Mass averaged SAR Setting:                                      */
    /************************************************************************/

	nXgSAR = 2;
	double XgSAR[] = 
    {/* Input mass you need compute(g/cm3) */
        1.0,
        10.0
    };

    /***END OF BLOCK 8.2*****************************************************/

/***END OF BLOCK 8***********************************************************/
#endif

/* End project setting */

/* Processing */
    strcpy(path_proj, set_path_proj);
    strcpy(proj_name, set_proj_name);
    strcpy(path_log, set_path_log);
    strcpy(path_save, set_path_save);

    strcpy(pathSRC, set_pathSRC);

    strcpy(path_data, set_path_data);
    strcpy(model_name, set_model_name);
    strcpy(media_name, set_media_name);

    fp_save_field_file = (field_file*) malloc (save_plane_amount * sizeof(field_file));/* Defaule setting: save plane in the plane where have feed point */
	for (i = 0; i < save_plane_amount; i++)
	{
		fp_save_field_file[i].sp.start=save_plane_parameters[i][0];
		fp_save_field_file[i].sp.end=save_plane_parameters[i][1];
		fp_save_field_file[i].sp.step=save_plane_parameters[i][2];
		fp_save_field_file[i].sp.plane_no=save_plane_parameters[i][3];
		fp_save_field_file[i].sp.slice=save_plane_parameters[i][4];
	}

#ifdef _SAR
    if (nXgSAR)
        save_localSAR_amount = _spaceZ;

   	pSAR = (localSAR*) calloc (save_localSAR_amount, sizeof(localSAR));
    if (nXgSAR == 0)
    {
	    for (i = 0; i<save_localSAR_amount; ++i)
	    {
	    	pSAR[i].start = save_sar_parameters[i][0];
	    	pSAR[i].end = save_sar_parameters[i][1];
	    	pSAR[i].plane_no = 1;
	    	pSAR[i].slice = save_sar_parameters[i][2];
	    }
    }
    else
    {
		for (i = 0; i<save_localSAR_amount; ++i)
        {
            pSAR[i].start = nMax-(int)(1/freq/dt)+1;
	    	pSAR[i].end = nMax;
	    	pSAR[i].plane_no = 1;
			pSAR[i].slice = i+1;
        }
    }
#endif
/* End processing */

/* Write a new project file */
if(myrank == 0)
{
    if (isNewProj == 0)
        return -1;
    char proj_path[MAX_SIZE_OF_PATH];
    strcpy(proj_path, path_proj);
    strcat(proj_path, proj_name);
    strcat(proj_path, ".txt");
    printf("Building new project --> %s.\n", proj_path);
	fflush(stdout);
	if (fp_proj = fopen(proj_path, "w+"))
	{
        fprintf(fp_proj, "<Project>\n");
		fprintf(fp_proj, "nMax=%d\n", nMax);
        fprintf(fp_proj, "dt=%e\n", dt);

        fprintf(fp_proj, "<Mesh>\n");
		fprintf(fp_proj, "_spaceX=%d,_spaceY=%d,_spaceZ=%d\n", _spaceX, _spaceY, _spaceZ);
        fprintf(fp_proj, "dx=%e,dy=%e,dz=%e\n", dx, dy, dz);
        fprintf(fp_proj, "padding=%d\n", padding);

        fprintf(fp_proj, "<Absorbing boundary>\n");
		fprintf(fp_proj, "abcNo=%d\n", abcNo);
		fprintf(fp_proj, "thicknessOfPml=%d\n", thicknessOfPml);

        fprintf(fp_proj, "<Power source>\n");
        fprintf(fp_proj, "sourceType=%d\n", sourceType);
		fprintf(fp_proj, "_isource=%d,_jsource=%d,_ksource=%d\n", _isource, _jsource, _ksource);
		fprintf(fp_proj, "port=%c\n", port);

		fprintf(fp_proj, "waveForm=%d\n", waveForm);
		fprintf(fp_proj, "amp=%lf\n", amp);

        fprintf(fp_proj, "freq=%lf\n", freq);
        fprintf(fp_proj, "t0=%d\n", t0);
        fprintf(fp_proj, "pulse_width=%d\n", pulse_width);
        fprintf(fp_proj, "pathSAR=%s\n", pathSRC);

        fprintf(fp_proj, "<Model>\n");
        fprintf(fp_proj, "-<Improt model>\n");
        fprintf(fp_proj, "path_data=%s\n", path_data);
        fprintf(fp_proj, "model_name=%s\n", model_name);
        fprintf(fp_proj, "media_name=%s\n", media_name);
		fprintf(fp_proj, "mediaNum=%d\n", mediaNum);
        fprintf(fp_proj, "-<Build object>\n");
        fprintf(fp_proj, "object_num=%d\n", object_num);
        for (i = 0; i<object_num; ++i)
        {
            fprintf(fp_proj, "%d,%d,%d,%d,%d,%d,%d,%lf,%lf,%lf\n", 3,
                    object_data[i][0], object_data[i][1], object_data[i][2], 
                    object_data[i][3], object_data[i][4], object_data[i][5], 
                    new_media[object_data[i][6]][0], 
                    new_media[object_data[i][6]][1], 
                    new_media[object_data[i][6]][2]);
        }
        fprintf(fp_proj, "-<Antenna>\n");
        fprintf(fp_proj, "antenna_amount=%d\n", antenna_amount);
        for (i = 0; i<antenna_amount; ++i)
        {
            fprintf(fp_proj, "antenna_direction=%d\n", antenna_direction);
            fprintf(fp_proj, "antenna_feed_x=%d,antenna_feed_y=%d,antenna_feed_z=%d\n", antenna_feed_x, antenna_feed_y, antenna_feed_z);
            fprintf(fp_proj, "antenna_impedance=%f\n", antenna_impedance);
            fprintf(fp_proj, "antenna_length_high=%d,antenna_length_low=%d\n", antenna_length_high, antenna_length_low);
        }

		fprintf(fp_proj, "<Field Save>\n");
		fprintf(fp_proj, "save_plane_amount=%d\n", save_plane_amount);
		for (i = 0; i < save_plane_amount; i++)
		{
			fprintf(fp_proj, "saveStart=%d,", fp_save_field_file[i].sp.start);
			fprintf(fp_proj, "saveEnd=%d,", fp_save_field_file[i].sp.end);
			fprintf(fp_proj, "saveStep=%d,", fp_save_field_file[i].sp.step);
			fprintf(fp_proj, "savePlaneNo=%d,", fp_save_field_file[i].sp.plane_no);
			fprintf(fp_proj, "slice=%d\n", fp_save_field_file[i].sp.slice);
		}
#ifdef _SAR
		fprintf(fp_proj, "<SAR>\n");

		fprintf(fp_proj, "-<Mass Averaged SAR>\n");
		fprintf(fp_proj, "Amount of Mass Averaged SAR %d\n", nXgSAR);
		for (i = 0; i<nXgSAR; ++i)
			fprintf(fp_proj, "Mass Averaged SAR %lf\n", XgSAR[i]);

		fprintf(fp_proj, "-<LocalSAR>\n");
		fprintf(fp_proj, "save_localSAR_amount=%d\n", save_localSAR_amount);
		for (i = 0; i < save_localSAR_amount; ++i)
		{
			fprintf(fp_proj, "saveLocalSARStart=%d,", pSAR[i].start);
			fprintf(fp_proj, "saveLocalSAREnd=%d,", pSAR[i].end);
			fprintf(fp_proj, "saveLocalSARPlaneNo=%d,", pSAR[i].plane_no);
			fprintf(fp_proj, "LocalSARslice=%d\n", pSAR[i].slice);
		}
#endif

		fclose(fp_proj);
		return SUCCESS;
	}
	else
		return FILE_ERROR;
}
else/* myrank != 0 */
    return SUCCESS;
}
