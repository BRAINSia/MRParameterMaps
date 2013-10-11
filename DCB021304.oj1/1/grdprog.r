;****************************************************************
;
; $Source: /pv/CvsTree/pv/gen/src/prg/parx_src/imnd/gp.tomo/gefi_tomo.r,v $
;
; Copyright (c) 1996
; Bruker Medizintechnik GmbH
; D-76275 Ettlingen, Germany
;
; All Rights Reserved
;
;
; $Locker:  $
; $State: Exp $
; $Log: gefi_tomo.r,v $
; Revision 3.6  1998/12/14 15:36:26  dwe
; CVS
;
; Revision 3.3.2.6  1998/12/11 13:16:32  dwe
;        DPH 111298: inversion module in tomo sequences
;        DWE 111298: spec package-  new diffusion codes- onepulse modules- reco_set
;                    in set_reco().
;
; Revision 3.3.2.5  1998/10/23 13:46:14  dwe
; *** empty log message ***
;
; Revision 3.5  1998/07/17 11:58:10  dwe
; *** empty log message ***
;
; Revision 3.4  1998/04/29 08:29:46  dwe
;
;
;****************************************************************

;written by DWE/HL 6.7.95
;read | phase | slice
loop NR;			Evolution-Loop
{
  loop l17 <3d>;	3D-Loop
  {
  loop ACQ_size[1] <2d>
  {
    begin dummy_projection
    loop NA;		Averaging-Loop
    {
	loop l26;		Movie-Loop
	{
	loop NSLICES <slice> <IMND_inv_matrix>
	{

if(IMND_suppression)
{
{     (-25,,no_scale)       |   (-25,,no_scale)      |  (-25,,no_scale)  }
{     (0)                   |   (0)                  |   (0)      }
{     (25,,no_scale)        |   (0)                  |   (0)      }
{     (0)                   |   (0)                  |   (0)      }
{     (0)                   |   (25,,no_scale)       |   (0)      }
{     (0)                   |   (0)                  |   (0)      }
{     (0)                   |   (0)                  |   (25,,no_scale)      }
{     (0)                   |   (0)                  |   (0)      }
}
if(IMND_fatflag)
{
{       (00)    |   (00)        |       (100)   }
{       (00)    |   (00)        |       (0)     }
}
if(IMND_invflag)
{
	loop l12 <IMND_InvRecov_grad>
	{
		{	(00)	|   (00)	|	(0),,,IMND_InvRecov_grad(100,,no_scale)}
	}
		{	(00)	|   (00)	|	(00)	}

}
if(IMND_FovSat_flag)
{
    loop l13 <IMND_FovSat_grad>
     {
{     (0),,,IMND_FovSat_grad(100,,direct_scale)       | (0)        | (0)  }
{     (0)       | (0)        | (0)+IMND_FovSat_spoiler(100,,no_scale)      }
     }
{     (0)       | (0)        | (0)      }
}
  if(IMND_InflowSat_flag)
  {
    loop L[14] <IMND_InflowSat_grad>
    {
	{     (0) |  (0)  | (0),,,IMND_InflowSat_grad(100,,no_scale)   }
	{     (0)       | (0)        | (0)+IMND_InflowSat_spoiler(100,,no_scale)      }
    }
   	{     (0)       | (0)        | (0)      }
  }

{     (t7,,no_scale)	| (0)	     | (0)	}
{     (0)	| (0)	     | (t0,,no_scale)	}
if(IMND_TE_long)
{
{     (0)		| (0)	     			| (0)	}
{     (0)		| (0)        			| (t1,,no_scale)}
{     (0)		| (0)	     			| (t1,,no_scale)}
{     (0)       	| (0)        			| (0)      }
{     (t4,,no_scale)	| (0),r2d(t2,,no_scale)	     	| (0),,r3d(t9,,no_scale)}
{     (t4,,no_scale)	| (0),r2d(t2,,no_scale)	     	| (0),,r3d(t9,,no_scale)}
{     (0)		| (0),r2d(t2,,no_scale)	     	| (0),,r3d(t9,,no_scale)}
}
if(!IMND_TE_long)
{
{     (t4,,no_scale)    | (0),r2d(t2,,no_scale)      	| (0),,r3d(t9,,no_scale)}
{     (t4,,no_scale)    | (0),r2d(t2,,no_scale)	     	| (t1,,no_scale),,r3d(t9,,no_scale)}
{     (0)       	| (0),r2d(t2,,no_scale)      	| (t1,,no_scale) ,,r3d(t9,,no_scale)}
}
{     (t5,,no_scale)    | (0)                           | (0)   }
if(IMND_grad_refo == Yes)
{
{     (t5,,no_scale)    | (0),r2d(t3,,no_scale)         | (0),,r3d(t8,,no_scale)}
}

{     (0)               | (0)                           | (t10,,no_scale)}

{     (0)		| (0)	     		     	| (0)	}


	};	Slice-Loop
	};	Movie-Loop
    };	Averaging-Loop
    end dummy_projection
};		2D-Loop
};		3D-Loop
};		Evolution-Loop
