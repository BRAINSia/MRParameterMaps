;msme_vtr2.r

loop NR
{
loop ACQ_size[1]/ACQ_rare_factor <2d> <rare>
{
  begin dummy_projection
  loop NA
  {
    loop NSLICES <slice>
    {
       { (80) | (0) |  (0)  }
       {  (0)  | (0) |  (0)  }
       {  (0)  | (0) | (80) }
       {  (0)  | (0) |  (0)  }
       { (t0)  | (0) | (t3)  }
       {  (0)  | (0) |  (0)  }

      if (ACQ_rare_factor ==1)
      {
      loop ACQ_rare_factor <rare>
      {
      loop NECHOES/ACQ_rare_factor
      {
       {  (0)  | (0) | (80) }
       {  (0)  | (0) |  (0)  }
     if(!IMND_hft_mode)
       {
       { (t1)  | (0),r2d(t6) | (t4)  }
       }
     if(IMND_hft_mode)
       {
       { (t1)  | (0),IMND_hft_grad(t6) | (t4)  }
       }
       {  (0)  | (0) |  (0)  }
       { (80) | (0) |  (0)  }
       {  (0)  | (0) |  (0)  }
     if(!IMND_hft_mode)
       {
       { (t2)  | (0),r2d(t7) | (t5)  }
       }
     if(IMND_hft_mode)
       {
       { (t2)  | (0),IMND_hft_grad(t6) | (t5) }
       }
       {  (0)  | (0) |  (0)  }
      }
      }
      }

      if (ACQ_rare_factor !=1)
      {
      loop NECHOES <rare>
      {
       {  (0)  | (0) | (80) }
       {  (0)  | (0) |  (0)  }
     if(!IMND_hft_mode)
       {
       { (t1)  | (0),r2d(t6) | (t4)  }
       }
     if(IMND_hft_mode)
       {
       { (t1)  | (0),IMND_hft_grad(t6) | (t4)  }
       }
       {  (0)  | (0) |  (0)  }
       { (80) | (0) |  (0)  }
       {  (0)  | (0) |  (0)  }
     if(!IMND_hft_mode)
       {
       { (t2)  | (0),r2d(t7) | (t5)  }
       }
     if(IMND_hft_mode)
       {
       { (t2)  | (0),IMND_hft_grad(t6) | (t5) }
       }
       {  (0)  | (0) |  (0)  }
      }
      }

      if (ACQ_flipback)
      {
       {  (0)  | (0) | (80) }
       {  (0)  | (0) |  (0)  }
       { (t0)  | (0) | (t3)  }
       {  (0)  | (0) |  (0)  }
       {  (0)  | (0) | (80) }
       {  (0)  | (0) |  (0)  }
      }
       {  (0)  | (0) |  (0)  }
  } } 
  end dummy_projection
}
}
