;****************************************************************
;
; $Source: /pv/CvsTree/pv/gen/src/prg/parx_src/imnd/gp.tomo/msme_mod.r,v $
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
; $Log: msme_mod.r,v $
; Revision 3.1  1997/09/30 14:00:07  gh
; Update from methods group
;
; Revision 1.1  1996/10/31  10:57:51  appadmin
; Autocheckin during RCSinit
;
;
;
;****************************************************************

;msme_mod.r
;Thu Dec 22 16:15:36 GMT 1994 KUS,KRE limit phase grad. duration

loop ACQ_size[1]/ACQ_rare_factor <2d> <rare>
{
  begin dummy_projection
  loop NA
  {
    loop NSLICES
    {
      {(t8)	|	(0)		|	(0)	}
      {(0)	|	(0)		|	(100)	}
      {(t0)	|	(0)		|	(t3)	}
      {(0)	|	(0)		|	(0)	};new
      loop NECHOES <rare>
      {
        {(0)	|	(0)		|	(100)	}
        {(0)	|	(0)		|	(0)	};new
        {(t1)	|	(0),r2d(t6)	|	(t4)	}
        {(t8)	|	(0)		|	(0)	}
        {(t2)	|	(0),r2d(t7)	|	(t5)	}
        {(0)	|	(0)		|	(0)	};new
      }
      if (ACQ_flipback)
      {
        {(0)	|	(0)		|	(100)	}
        {(t0)	|	(0)		|	(t3)	}
        {(0)	|	(0)		|	(0)	};new
        {(0)	|	(0)		|	(100)	}
      }
      {(0)	|	(0)		|	(100)	}
      {(0)	|	(0)		|	(0)	}
  } } 
  end dummy_projection
}

