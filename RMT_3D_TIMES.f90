MODULE RMT_3D_TIMES
PUBLIC :: TimeS
CONTAINS
SUBROUTINE TimeS(The_Time)
integer The_Time
integer,dimension(8) :: vals
call DATE_AND_TIME(VALUES=vals)
The_Time  = vals(5)*3600 + vals(6)*60 + vals(7)
END SUBROUTINE TimeS
END MODULE RMT_3D_TIMES