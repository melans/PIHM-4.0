#include "IO.hpp"

void mkdir_p( char *dir, int mode) {
    /* MODIFIED FROM http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html */
    char tmp[MAXLEN];
    char *p = NULL;
    size_t len;
    
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if(tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            mkdir(tmp, mode);
            *p = '/';
        }
    mkdir(tmp, mode);
}
void FileIn::saveProject(){
    FILE *fp;
    char fn[MAXLEN];
    sprintf(fn, "%s/%s.prj", outpath, projectname);
    fp = fopen(fn,"w");
    CheckFile(fp, fn);
    fprintf(fp, "PRJ \t %s\n", projectname);
    fprintf(fp, "INPATH \t %s\n", inpath);
    fprintf(fp, "OUTPATH \t %s\n", outpath);
//    if(ilog){
//        fprintf(fp, "LOG \t %s\n", logfile);
//    }
    fprintf(fp, "MESH \t %s\n", file_mesh);
    fprintf(fp, "ATT \t %s\n", file_att);
    fprintf(fp, "INIT \t %s\n", file_init);
    fprintf(fp, "CALIB \t %s\n", file_calib);
    fprintf(fp, "RIV \t %s\n", file_riv);
    fprintf(fp, "RIV \t %s\n", file_rivchn);
    fprintf(fp, "LC \t %s\n", file_lc);
    fprintf(fp, "SOIL \t %s\n", file_soil);
    fprintf(fp, "GEOL \t %s\n", file_geol);
    fprintf(fp, "FORC \t %s\n", file_forc);
    fprintf(fp, "LAI \t %s\n", file_lai);
    fprintf(fp, "MF \t %s\n", file_mf);
    fprintf(fp, "IBC \t %s\n", file_bc);
    fprintf(fp, "PARA \t %s\n", file_para);
    fprintf(fp, "LCM \t %s\n", file_lcm);
    fclose(fp);
}

void FileIn:: setInFilePath(char * indir, char *  pjrname, int nthred){
    numthreads = nthred;
    setInFilePath(indir, pjrname);
}
void FileIn:: setInFilePath(char * indir, char *  pjrname){
    /*Model io*/
    sprintf(inpath, "%s", indir);
    sprintf(outpath, "output/%s.out", pjrname);
    sprintf(logfile, "%s.log", outpath);
    sprintf(projectname, "%s",  pjrname);
    
    /*Spatial Data*/
    sprintf(file_mesh, "%s/%s.%s", inpath, projectname, "sp.mesh");
    sprintf(file_att, "%s/%s.%s", inpath, projectname, "sp.att");
    sprintf(file_riv, "%s/%s.%s", inpath, projectname, "sp.riv");
    sprintf(file_rivchn, "%s/%s.%s", inpath, projectname, "sp.rivchn");
    sprintf(file_lake, "%s/%s.%s", inpath, projectname, "sp.lake");
    
    /* physical parameters */
    sprintf(file_lc, "%s/%s.%s", inpath, projectname, "para.lc");
    sprintf(file_soil, "%s/%s.%s", inpath, projectname, "para.soil");
    sprintf(file_geol, "%s/%s.%s", inpath, projectname, "para.geol");
    
    /* model configuration */
    sprintf(file_para, "%s/%s.%s", inpath, projectname, "cfg.para");
    sprintf(file_calib, "%s/%s.%s", inpath, projectname, "cfg.calib");
    sprintf(file_init, "%s/%s.%s", inpath, projectname, "cfg.ic");
    
    /* Time-series data */
    sprintf(file_forc, "%s/%s.%s", inpath, projectname, "tsd.forc");
    sprintf(file_lai, "%s/%s.%s", inpath, projectname, "tsd.lai");
    sprintf(file_mf, "%s/%s.%s", inpath, projectname, "tsd.mf");
    sprintf(file_rl, "%s/%s.%s", inpath, projectname, "tsd.rl");
    sprintf(file_bc, "%s/%s.%s", inpath, projectname, "tsd.bc");
    sprintf(file_lcm, "%s/%s.%s", inpath, projectname, "tsd.lcm");
    sprintf(file_obs, "%s/%s.%s", inpath, projectname, "tsd.obs");
}

FileOut::FileOut(){
    setsuffix("");
}
void FileOut::setsuffix(const char *s){
    strcpy(suffix, s);
}
void FileOut::setOutFilePath(char *outdir, char *  prjname){
    sprintf(outpath, "%s", outdir);
    sprintf(projectname, "%s", prjname);
    updateFilePath();
}
void FileOut::updateFilePath(){
    createDir();
    //    printf("\nOutput path is \"%s\"\n", outpath);
    /****************************************************
     format： Prjname+Surffix.Component+Type+Value.Extension
     i.g. pj0.eleysurf.dat/ prject001.rivqdown.dat
     pj - project name. len = 1~n
     0  - suffix. len = 1~n
     ele- componencts, len = 3
            riv = River
            ele = Element
            lak = Lake
     y  - Type,
            y - State Variable[m],
            v- specific flux [m3/m2/day],
            q - flux variable[m3/day]; len =1
     surf- variable, len = 1~n;
     dat- File extension, len =3
            dat = bindary,
            csv = ASCII spreadsheet.
     *****************************************************/
    //rivers
    sprintf(riv_Q_down, "%s/%s%s.rivqdown", outpath, projectname, suffix);
    sprintf(riv_Q_up, "%s/%s%s.rivqup", outpath, projectname, suffix);
    sprintf(riv_Q_surf, "%s/%s%s.rivqsurf", outpath, projectname, suffix);
    sprintf(riv_Q_sub, "%s/%s%s.rivqsub", outpath, projectname, suffix);
    sprintf(riv_y_stage, "%s/%s%s.rivystage", outpath, projectname, suffix);
    //cells
    sprintf(ele_y_snow, "%s/%s%s.eleysnow", outpath, projectname, suffix);
    sprintf(ele_y_ic, "%s/%s%s.eleyic", outpath, projectname, suffix);
    sprintf(ele_y_surf, "%s/%s%s.eleysurf", outpath, projectname, suffix);
    sprintf(ele_y_unsat, "%s/%s%s.eleyunsat", outpath, projectname, suffix);
    sprintf(ele_y_wetfrount, "%s/%s%s.eleywf", outpath, projectname, suffix);
    sprintf(ele_y_gw, "%s/%s%s.eleygw", outpath, projectname, suffix);
    //cell-fluxes
    sprintf(ele_q_ET[0], "%s/%s%s.elevetic", outpath, projectname, suffix);
    sprintf(ele_q_ET[1], "%s/%s%s.elevettr", outpath, projectname, suffix);
    sprintf(ele_q_ET[2], "%s/%s%s.elevetev", outpath, projectname, suffix);
    
    //    sprintf(ele_Q_surf, "%s/%s%s.eleqsurf", outpath, projname);
    //    sprintf(ele_Q_sub, "%s/%s%s.eleqsub", outpath, projname);
    sprintf(ele_Q_subTot, "%s/%s%s.eleqsub", outpath, projectname, suffix);
    sprintf(ele_Q_surfTot, "%s/%s%s.eleqsurf", outpath, projectname, suffix);
    
    sprintf(ele_q_ETP, "%s/%s%s.elevetp", outpath, projectname, suffix);
    sprintf(ele_q_prcp, "%s/%s%s.elevprcp", outpath, projectname, suffix);
    sprintf(ele_q_netprcp, "%s/%s%s.elevnetprcp", outpath, projectname, suffix);
    sprintf(ele_q_infil, "%s/%s%s.elevinfil", outpath, projectname, suffix);
    sprintf(ele_q_rech, "%s/%s%s.elevrech", outpath, projectname, suffix);
    //cell_wb
    sprintf(ewb_q_in, "%s/%s%s.ewbqin", outpath, projectname, suffix);
    sprintf(ewb_q_out, "%s/%s%s.ewbqout", outpath, projectname, suffix);
    sprintf(ewb_q_io, "%s/%s%s.ewbqio", outpath, projectname, suffix);
    sprintf(ewb_dh, "%s/%s%s.ewbydh", outpath, projectname, suffix);
    
    sprintf(lake_Q_riv, "%s/%s%s.lakqriv", outpath, projectname, suffix);
    sprintf(lake_Q_surf, "%s/%s%s.lakqsurf", outpath, projectname, suffix);
    sprintf(lake_Q_sub, "%s/%s%s.lakqsub", outpath, projectname, suffix);
    sprintf(lake_y_stage, "%s/%s%s.lakystage", outpath, projectname, suffix);
    sprintf(lake_q_evap, "%s/%s%s.lakvevap", outpath, projectname, suffix);
    sprintf(lake_q_prcp, "%s/%s%s.lakvprcp", outpath, projectname, suffix);
    
    sprintf(Init_update, "%s/%s%s.update.ic", outpath, projectname, suffix);
    sprintf(Init_bak, "%s/%s%s.cfg.ic.bak", outpath, projectname, suffix);
    
    sprintf(Calib_bak, "%s/%s%s.cfg.calib.bak", outpath, projectname, suffix);
    
    sprintf(floodout, "%s/%s%s.flood.csv", outpath, projectname, suffix);
    sprintf(obs_sim, "%s/%s%s.ovs.csv", outpath, projectname, suffix);
}
void FileOut::createDir(){
    mkdir_p(outpath, 0777);
}
void FileIn::setOutpath(const char *fn){
    strcpy(outpath, fn);
}
void FileIn::setCalibFile(const char *fn){
    strcpy(file_calib, fn);
}
void FileIn::readProject(const char *fn){
    FILE *fp;
    char *p;
    int i, j;
    int ip;
    int num = 17;
    char keywords[17][MAXLEN] = {"PRJ", "inpath", "outpath", "mesh", "att", "riv", "para", "calib", "lc",
        "soil", "geol", "forc", "lai", "ibc", "init","lcm", "mf"};
    char key[MAXLEN];
    char value[MAXLEN];
    char empty[MAXLEN];
    char str[MAXLEN];
    char prjname[MAXLEN];
    char inpath[MAXLEN];
    char outpath[MAXLEN];
    fp = fopen(fn,"r");
    CheckFile(fp, fn);
    strcpy(file_lcm, "");
    ip = -1;
    fgets(str,MAXLEN, fp);
    for(i = 0; i < num && !feof(fp) ; i++){
        
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0')
        {
            continue;
        }
        
        sscanf(str, "%s %s",key, value);
        for(j = 0; j < num; j++) {
            if(strcasecmp(key, keywords[j])==0){
                ip = j;
                break;
            }
        }
        switch (ip) {
            case 0    :
                p = projectname;
                break;
            case 1    :
                p = inpath;
                break;
            case 2    :
                p = outpath;
                break;
            case 3    :
                p = file_mesh;
                break;
            case 4    :
                p = file_att;
                break;
            case 5    :
                p = file_riv;
                sprintf(file_rivchn, "%schn", file_riv);
                break;
            case 6    :
                p = file_para;
                break;
            case 7    :
                p = file_calib;
                break;
            case 8    :
                p = file_lc;
                break;
            case 9    :
                p = file_soil;
                break;
            case 10   :
                p = file_geol;
                break;
            case 11   :
                p = file_forc;
                break;
            case 12   :
                p = file_lai;
                break;
            case 13   :
                p = file_bc;
                break;
            case 14   :
                p = file_init;
                break;
            case 15   :
                p = file_lcm;
                break;
            case 16   :
                p = file_mf;
                break;
            default:
                p = empty;
                break;
        }
        if(ip ==  0){
            strcpy(prjname, projectname);
            sprintf(inpath, "input/%s", prjname);
            sprintf(outpath, "output/%s.out", prjname);
            setInFilePath(inpath, prjname);
            strcpy(outpath, outpath);
        }
        strcpy(p, value);
        fgets(str,MAXLEN, fp);
    }
    
    fclose(fp);
}
void FileOut::setOutpath(const char *fn){
    strcpy(outpath, fn);
    updateFilePath();
}
void FileOut::copy(FileOut *p){
    strcpy(   outpath, p->outpath );
    strcpy(   projectname, p->projectname );
    strcpy(   riv_Q_down, p->riv_Q_down );
    strcpy(   riv_Q_surf, p->riv_Q_surf );
    strcpy(   riv_Q_sub, p->riv_Q_sub );
    strcpy(   riv_y_stage, p->riv_y_stage );
    strcpy(   ele_y_snow, p->ele_y_snow );
    strcpy(   ele_y_ic, p-> ele_y_ic);
    strcpy(   ele_y_surf, p->ele_y_surf );
    strcpy(   ele_y_unsat, p->ele_y_unsat );
    strcpy(   ele_y_gw, p->ele_y_gw );
    strcpy(   ele_y_wetfrount, p->ele_y_wetfrount );
    for(int i = 0; i < 3; i++){
        strcpy(   ele_q_ET[i], p->ele_q_ET[i] );
    }
    strcpy(   ele_q_ETP, p->ele_q_ETP );
    strcpy(   ele_q_prcp, p->ele_q_prcp );
    strcpy(   ele_q_netprcp, p->ele_q_netprcp );
    strcpy(   ele_Q_surfTot, p->ele_Q_surfTot );
    strcpy(   ele_Q_subTot, p->ele_Q_subTot );
    strcpy(   ele_q_infil, p->ele_q_infil );
    strcpy(   ele_q_rech, p->ele_q_rech );
    
    strcpy(   ewb_q_in, p->ewb_q_in );
    strcpy(   ewb_q_out, p->ewb_q_out );
    strcpy(   ewb_q_io, p->ewb_q_io );
    strcpy(   ewb_dh, p->ewb_dh );
    
    strcpy(   lake_Q_riv, p->lake_Q_riv );
    strcpy(   lake_Q_surf, p->lake_Q_surf );
    strcpy(   lake_Q_sub, p->lake_Q_sub );
    strcpy(   lake_y_stage, p->lake_y_stage );
    strcpy(   lake_q_evap, p->lake_q_evap );
    strcpy(   lake_q_prcp, p->lake_q_prcp );
    
    strcpy(   Init_update, p->Init_update );
    strcpy(   Init_bak, p->Init_bak );
    
    strcpy(   floodout, p->floodout );
    strcpy(   obs_sim, p->obs_sim );

}

