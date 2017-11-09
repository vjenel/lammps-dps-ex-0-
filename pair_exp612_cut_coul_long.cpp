/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include "pair_lj_cut_coul_long.h" // In original code 
#include "pair_exp612_cut_coul_long.h"  //added
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairExp612CutCoulLong::PairExp612CutCoulLong(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  respa_enable = 1;
  writedata = 1;
  ftable = NULL;
  qdist = 0.0;
}

/* ---------------------------------------------------------------------- */

PairExp612CutCoulLong::~PairExp612CutCoulLong()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
/* In original code of lj/cut/coul/long:
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
*/
    memory->destroy(AA); //added
    memory->destroy(bb); //added
    memory->destroy(CC6); //added
    memory->destroy(DD12); //added
    memory->destroy(offset); //added
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairExp612CutCoulLong::compute(int eflag, int vflag)
{
  int i,ii,j,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,r_bb,rexp,r2inv,r12inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq,b_b;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;



  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];


      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
//if (i==0||j==0){
//printf(" ++++++++ before mask i j %d %d factor_lj factor_coul %f %f sbmask(j) %d \n", i+1,j+1, factor_lj, factor_coul, sbmask(j));
//}

      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          b_b = bb[itype][jtype];
          r_bb = r*b_b;
          rexp = exp(-r_bb);
          r6inv = r2inv*r2inv*r2inv;
          r12inv = r6inv*r6inv;
          //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); // this was in original code of pait lj/cut/coul/long
          forcelj =   (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ; //added 
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            //evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);// in original code of pait lj/cut/coul/long
            evdwl = AA[itype][jtype]*rexp - r6inv*CC6[itype][jtype] + r12inv*DD12[itype][jtype]- offset[itype][jtype]; // added
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairExp612CutCoulLong::compute_inner()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double r,b_b,r_bb,rsq,rexp,r2inv,r6inv,r12inv,forcecoul,forcelj,factor_coul,factor_lj;
  double rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = listinner->inum;
  ilist = listinner->ilist;
  numneigh = listinner->numneigh;
  firstneigh = listinner->firstneigh;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms
//printf("  I am in rirLJCutCoulLong::compute_inner inner \n", i,j);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq) {
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        jtype = type[j];
        if (rsq < cut_ljsq[itype][jtype]) {
          r=sqrt(rsq);
          b_b = bb[itype][jtype];
          r_bb = r*b_b;
          rexp = exp(-r_bb);
          r6inv = r2inv*r2inv*r2inv;
          r12inv=r6inv*r6inv; 
          // forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); // original code
          forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ;
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairExp612CutCoulLong::compute_middle()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double r,rsq,r12inv,b_b,r_bb,r2inv,r6inv,rexp,forcecoul,forcelj,factor_coul,factor_lj;
  double rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = listmiddle->inum;
  ilist = listmiddle->ilist;
  numneigh = listmiddle->numneigh;
  firstneigh = listmiddle->firstneigh;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms
//printf("  I am in PairLJCutCoulLong::compute_middle \n", i,j);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      r2inv = 1.0/rsq;
      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        jtype = type[j];
        if (rsq < cut_ljsq[itype][jtype]) {
          r=sqrt(rsq);
          b_b = bb[itype][jtype];
          r_bb = r*b_b;
          rexp = exp(-r_bb);
          r6inv = r2inv*r2inv*r2inv;
          r12inv = r6inv*r6inv;
          // forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);  // original code
          forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ;
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        if (rsq < cut_in_on_sq) {
          rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
          fpair *= rsw*rsw*(3.0 - 2.0*rsw);
        }
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairExp612CutCoulLong::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rexp,r_bb,r2inv,r12inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj, b_b;
  double grij,expm2,prefactor,t,erfc;
  double rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = listouter->inum;
  ilist = listouter->ilist;
  numneigh = listouter->numneigh;
  firstneigh = listouter->firstneigh;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  // loop over neighbors of my atoms
//printf("  I am in PairLJCutCoulLong::compute_outer \n", i,j);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
/// I calculate these here because I dont really plan to use tables hence it is more efficient to compuytre them here
        r = sqrt(rsq);  
        r6inv = r2inv*r2inv*r2inv;
        r12inv = r6inv*r6inv;  //added
          b_b = bb[itype][jtype]; // added
          r_bb = r*b_b;
          rexp = exp(-r_bb);     // added

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - 1.0);
            if (rsq > cut_in_off_sq) {
              if (rsq < cut_in_on_sq) {
                rsw = (r - cut_in_off)/cut_in_diff;
                forcecoul += prefactor*rsw*rsw*(3.0 - 2.0*rsw);
                if (factor_coul < 1.0)
                  forcecoul -=
                    (1.0-factor_coul)*prefactor*rsw*rsw*(3.0 - 2.0*rsw);
              } else {
                forcecoul += prefactor;
                if (factor_coul < 1.0)
                  forcecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype] && rsq > cut_in_off_sq) {
          // forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); // in the original code of pair lj/cut/coul/long 
          forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ;
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            forcelj *= rsw*rsw*(3.0 - 2.0*rsw);
          }
        } else forcelj = 0.0;

        fpair = (forcecoul + forcelj) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq) {
              ecoul = prefactor*erfc;
              if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
            } else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
              if (factor_coul < 1.0) {
                table = ptable[itable] + fraction*dptable[itable];
                prefactor = qtmp*q[j] * table;
                ecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            //evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -offset[itype][jtype];//In the original code
            evdwl = AA[itype][jtype]*rexp - r6inv*CC6[itype][jtype] + r12inv*DD12[itype][jtype]- offset[itype][jtype]; // added
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (vflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq) {
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
              if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
            } else {
              table = vtable[itable] + fraction*dvtable[itable];
              forcecoul = qtmp*q[j] * table;
              if (factor_coul < 1.0) {
                table = ptable[itable] + fraction*dptable[itable];
                prefactor = qtmp*q[j] * table;
                forcecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else forcecoul = 0.0;

          if (rsq <= cut_in_off_sq) {
             //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);// in the original code
             forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ;
          } else if (rsq <= cut_in_on_sq) {
             // forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); // in the original code
             forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ;
          }
          fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,ecoul,fpair,delx,dely,delz);
      } // if (rsq < cutsq[itype][jtype]) {
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");

/* These cpommented lines were in the original code of the pait lj/cut/coul/long
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
*/
  memory->create(AA,n+1,n+1,"pair:epsilon"); //added
  memory->create(bb,n+1,n+1,"pair:sigma");   // added
  memory->create(CC6,n+1,n+1,"pair:lj1");    // added
  memory->create(DD12,n+1,n+1,"pair:lj2");   // added
  memory->create(offset,n+1,n+1,"pair:offset");  //added
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::settings(int narg, char **arg)
{
 if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

/* The original code of PairLJCutCoulLong::coeff listed below was modified as in 
// PairExp612CutCoulLong::coeff
void PairLJCutCoulLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_lj_one = cut_lj_global;
  if (narg == 5) cut_lj_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
*/

void PairExp612CutCoulLong::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double AA_one = force->numeric(FLERR,arg[2]);
  double bb_one = force->numeric(FLERR,arg[3]);
  double CC6_one = force->numeric(FLERR,arg[4]);
  double DD12_one = force->numeric(FLERR,arg[5]);

  double cut_lj_one = cut_lj_global;
//  if (narg == 6) cut_lj_one = force->numeric(FLERR,arg[6]);
printf(" %d %d  %d  %lf %lf %lf %lf \n", ilo,jlo, narg, AA_one,bb_one,CC6_one,DD12_one);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      AA[i][j] = AA_one;
      bb[i][j] = bb_one;
      CC6[i][j] = CC6_one;
      DD12[i][j] = DD12_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style exp612/cut/coul/long requires atom attribute q");

  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this,instance_me);
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,cut_respa);
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */


/*
// The original code of PairLJCutCoulLong::init_one listed(commented) below
// was modified as in PairExp612CutCoulLong::init_one
double PairLJCutCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
  }

  // include TIP4P qdist in full cutoff, qdist = 0.0 if not TIP4P

  double cut = MAX(cut_lj[i][j],cut_coul+2.0*qdist);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut;
}

*/

double PairExp612CutCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
//    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
//                               sigma[i][i],sigma[j][j]);
//    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
  }

  // include TIP4P qdist in full cutoff, qdist = 0.0 if not TIP4P

  double cut = MAX(cut_lj[i][j],cut_coul+2.0*qdist);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double c = cut_lj[i][j];
    double c6=pow(1.0/c,6);
    double c12=c6*c6;
    double cexp = exp(-c*bb[i][j]);
    offset[i][j] = AA[i][j]*cexp - c6*CC6[i][j] + c12*DD12[i][j];
  } else offset[i][j] = 0.0;

     cut_ljsq[j][i] = cut_ljsq[i][j];
     AA[j][i] = AA[i][j];
     bb[j][i] = bb[i][j];
     CC6[j][i] = CC6[i][j];
     DD12[j][i] = DD12[i][j];
     offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double c=cut_lj[i][j], c3=c*c*c, c9=c3*c3, b_b_3=bb[i][j]*bb[i][j]*bb[i][j], b_b_4=b_b_3*bb[i][j];
    double bbc = bb[i][j]*c, bbc2 = bbc*bbc;
    double e_x_p = exp(-bbc);
    
    double t1 = AA[i][j]*(bbc2+2.0*bbc+2.0)*e_x_p;
    double t2 = AA[i][j] * e_x_p * bb[i][j] * (bbc*bbc2+3.0*bbc2+6.0*bbc+6.0);
    double all_01 = all[0]*all[1] * (2.0*MY_PI) ; 
//    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *  sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);   // in original code
//    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] * sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9); // in original code

    etail_ij = all_01*(t1/b_b_3 - CC6[i][j]/(3.0*c3) + DD12[i][j]/(9.0*c9));
    ptail_ij = 1.0/3.0*all_01*(t2/b_b_4 - 2.0*CC6[i][j]/(c3) + 4.0*DD12[i][j]/(3.0*c9));
  }

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
//        fwrite(&epsilon[i][j],sizeof(double),1,fp); // original code of lj/cut/coul/long
//        fwrite(&sigma[i][j],sizeof(double),1,fp);   // original code of lj/cut/coul/long
//        fwrite(&cut_lj[i][j],sizeof(double),1,fp);  // original code of lj/cut/coul/long

        fwrite(&AA[i][j],sizeof(double),1,fp);  //added
        fwrite(&bb[i][j],sizeof(double),1,fp);  //added
        fwrite(&CC6[i][j],sizeof(double),1,fp);  //added
        fwrite(&DD12[i][j],sizeof(double),1,fp);  //added
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);  //added
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
//          fread(&epsilon[i][j],sizeof(double),1,fp); // original code of lj/cut/coul/long
//          fread(&sigma[i][j],sizeof(double),1,fp);   // original code of lj/cut/coul/long
//          fread(&cut_lj[i][j],sizeof(double),1,fp);  // original code of lj/cut/coul/long

          fread(&AA[i][j],sizeof(double),1,fp); //added
          fread(&bb[i][j],sizeof(double),1,fp);  //added
          fread(&CC6[i][j],sizeof(double),1,fp);  //added
          fread(&DD12[i][j],sizeof(double),1,fp);  //added
          fread(&cut_lj[i][j],sizeof(double),1,fp); //added
        }
        MPI_Bcast(&AA[i][j],1,MPI_DOUBLE,0,world);  //added
        MPI_Bcast(&bb[i][j],1,MPI_DOUBLE,0,world);  //added
        MPI_Bcast(&CC6[i][j],1,MPI_DOUBLE,0,world);  //added
        MPI_Bcast(&DD12[i][j],1,MPI_DOUBLE,0,world);  //added
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);  //added
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
    fread(&ncoultablebits,sizeof(int),1,fp);
    fread(&tabinner,sizeof(double),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
//     fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]); // in the original code 
    fprintf(fp,"%d %g %g %g %g \n",i,AA[i][i],bb[i][i],CC6[i][i],DD12[i][i]);  // added
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairExp612CutCoulLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) 
//      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut_lj[i][j]); // in the original code
      fprintf(fp,"%d %d %g %g %g %g %g \n",i,j,AA[i][j],bb[i][j],CC6[i][j],DD12[i][j],cut_lj[i][j]); // added
}

/* ---------------------------------------------------------------------- */

double PairExp612CutCoulLong::single(int i, int j, int itype, int jtype,
                                 double rsq,
                                 double factor_coul, double factor_lj,
                                 double &fforce)
{
  double r2inv,r6inv,r,r_bb,r12inv,grij,expm2,t,erfc,prefactor,b_b,rexp;
  double fraction,table,forcecoul,forcelj,phicoul,philj;
  int itable;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      t = 1.0 / (1.0 + EWALD_P*grij);
      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
      prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    } else {
      union_int_float_t rsq_lookup_single;
      rsq_lookup_single.f = rsq;
      itable = rsq_lookup_single.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup_single.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction*dftable[itable];
      forcecoul = atom->q[i]*atom->q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction*dctable[itable];
        prefactor = atom->q[i]*atom->q[j] * table;
        forcecoul -= (1.0-factor_coul)*prefactor;
      }
    }
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
    r_bb = r*bb[itype][jtype];
    rexp = exp(-r_bb);
    r6inv = r2inv*r2inv*r2inv;
    r12inv = r6inv*r6inv;
    // forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]); // original code
    forcelj = (r_bb*AA[itype][jtype])*rexp - 6.0*r6inv*CC6[itype][jtype] + 12.0*r12inv*DD12[itype][jtype] ; // added
  } else forcelj = 0.0;

  fforce = (forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq)
      phicoul = prefactor*erfc;
    else {
      table = etable[itable] + fraction*detable[itable];
      phicoul = atom->q[i]*atom->q[j] * table;
    }
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }

  if (rsq < cut_ljsq[itype][jtype]) {
    // philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -offset[itype][jtype]; // original code
    philj = AA[itype][jtype]*rexp -  r6inv*CC6[itype][jtype] + r12inv*DD12[itype][jtype] - offset[itype][jtype]; // added
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairExp612CutCoulLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
//  dim = 2;   // original code
  //if (strcmp(str,"epsilon") == 0) return (void *) epsilon;  // original code
  //  if (strcmp(str,"sigma") == 0) return (void *) sigma;  // original code

  dim = 2;  ///added  (left like in the original code)
  if (strcmp(str,"AA") == 0) return (void *) AA; //added
  if (strcmp(str,"bb") == 0) return (void *) bb; //added
  if (strcmp(str,"CC6") == 0) return (void *) CC6;  //added
  if (strcmp(str,"DD12") == 0) return (void *) DD12; //added
  return NULL;
}
