# RTK_Construct
  ตัวอย่างการออกแบบระบบพิกัดความคลาดเคลื่อนต่ำสำหรับพื้นที่โครงการ “บ่อพลอย” ได้จากการพิจารณากำหนดจุดสุ่มทดสอบใช้หมุดควบคุม RID-GNSS-001 ถึง RID-GNSS-176 ทั้งหมดในโครงการ พร้อมจุดตัวแทนแนวกึ่งกลางของอุโมงค์ส่งน้ำ T1, T2 และ T3 กำหนดค่าออกแบบเบื้องต้นเป็น
Point	UTM_N	UTM_E	MSL
Tunnel-Inlet	1,608,603.000	518,123.000	160.0
Tunnel-Mid	1,608,118.000	529,367.000	158.5
Tunnel-Outlet	1,607,439.000	539,639.000	157.0

  ความสูงเหนือทะเลปานกลางของหมุดควบคุมได้จากงานเดินระดับดิฟเฟอเรนเชียลงานชั้น 2 โครงการเลือกใช้การฉายแผนที่ระบบพิกัดทรานเวอร์สเมอร์เคเตอร์ (Transverse Mercator :TM ) ให้ศูนย์กลางเมอริเดียนเป็น 99°30’ การคำนวณไปความสูงเหนือทรงรีใช้แบบจำลองยีออยด์ TGM-2017 เนื่องจากมีงานออกแบบและก่อสร้างอยูลึกใต้ผิวภูมิประเทศลงไป จึงมีการกำหนดจุดทดสอบ T1, T2 และ T3 ลึกใต้ผิวภูมิประเทศ ในขั้นต้นระดับ projection plane อยู่ในเกณฑ์ค่าเฉลี่ยของระดับทะเลปานกลางของจุดทดสอบ (offset PP=0 m) ได้ค่าเสกลแฟกเตอร์ร่วม CSF min /max : -12.3 / 24.7 ppm ซึ่งเกินเกิน +20 ppm ดังนั้นเพื่อให้ระนาบการฉายแผนที่อยู่ในเกณฑ์ค่าระดับของอุโมงค์ด้วยจึงเปลี่ยนระดับการฉายต่ำลง offset PP = -20 m พบว่าได้ผลเป็นที่น่าใจ  ได้ค่าเสกลแฟกเตอร์ร่วม CSF min /max -18.3 / +18.7 ppm
  
Offset_PP	Points MSL min/mean/max : +30.7 / +76.6 / +202.6 meter
E-W : 70,016 m. 		 N-S : 59,889 m.
+20 m	Designed Project Plane : hPP=+62.54 m. / MSL=+96.39 m. equiv. to k0 = 1.000010
Points CSF min/mean/max : -12.3 / +8.7 / +24.7 ppm
0 m	Designed Project Plane : hPP=+42.54 m. / MSL=+76.39 m. equiv. to k0 = 1.000007
Points CSF min/mean/max : -15.3 / +5.7 / +21.7 ppm
-20 m	Designed Project Plane : hPP=+22.54 m. / MSL=+56.39 m. equiv. to k0 = 1.000004
Points CSF min/mean/max : -18.3 / +2.7 / +18.7 ppm
-40 m	Designed Project Plane : hPP=+2.54 m. / MSL=+36.39 m. equiv. to k0 = 1.000000
Points CSF min/mean/max : -22.3 / -1.3 / +14.7 ppm

 
ทำการปรับค่า False Easting และ Fase Northing ให้พิกัดท้องถิ่นมีตัวเลขกระทัดรัดและกระชับ ไม่มีตัวเลขค่าพิกัดมีเครื่องหมายลบ ในมาตรฐานฉบับนี้จึงเสนอว่า ค่าพารามิเตอร์ระบบพิกัดท้องถิ่น LDP ของโครงการ “บ่อพลอย” (LDP-BOPLOI) คือ 
ค่าพารามิเตอร์ระบบพิกัดท้องถิ่น LDP ของโครงการ “บ่อพลอย” (LDP-BOPLOI)

+proj=tmerc +lat_0=0.0 +lon_0=99.5 +k_0=1.000004  +x_0=50000  +y_0=-1550000
+a=6378137.0 +b=6356752.314245179 +units=m +no_defs +type=crs
+proj=tmerc +lat_0=0.0 +lon_0=99°30′ +k_0=1.000004  +x_0=50000  +y_0=-1550000
+a=6378137.0 +b=6356752.314245179 +units=m +no_defs +type=crs

Designed Project Plane : hPP=+22.54 m. / MSL=+56.39 m. equiv. to k0 = 1.000004

|     | Point         |       lng |       lat |   MSL |   CSF_ppm |      LDP_E |      LDP_N |
|----:|:--------------|----------:|----------:|------:|----------:|-----------:|-----------:|
|   0 | Tunnel-Inlet  | 99.168228 | 14.550374 |   160 |       0.1 | 14,243.947 | 59,272.457 |
|   1 | Tunnel-Mid    | 99.272595 | 14.545890 |   158 |      -8.1 | 25,491.449 | 58,762.602 |
|   2 | Tunnel-Outlet | 99.367932 | 14.539621 |   157 |     -12.8 | 35,765.983 | 58,060.819 |
|   3 | RID-GNSS-001  | 99.168657 | 14.552301 |   185 |      -3.9 | 14,290.498 | 59,485.635 |
|   4 | RID-GNSS-002  | 99.168475 | 14.551096 |   188 |      -4.3 | 14,270.748 | 59,352.315 |
|   5 | RID-GNSS-003  | 99.360889 | 14.546048 |   193 |     -18.1 | 35,007.322 | 58,772.421 |
|   6 | RID-GNSS-004  | 99.359894 | 14.545350 |   194 |     -18.3 | 34,900.101 | 58,695.172 |
|   7 | RID-GNSS-005  | 99.375530 | 14.539136 |   161 |     -13.6 | 36,584.932 | 58,006.718 |
|   8 | RID-GNSS-006  | 99.374453 | 14.537768 |   157 |     -13.1 | 36,468.759 | 57,855.352 |
|   9 | RID-GNSS-007  | 99.420384 | 14.521585 |   134 |     -10.8 | 41,418.474 | 56,062.594 |
|  10 | RID-GNSS-008  | 99.420039 | 14.520434 |   132 |     -10.4 | 41,381.274 | 55,935.318 |
|  11 | RID-GNSS-009  | 99.454297 | 14.491703 |   114 |      -8.2 | 45,073.232 | 52,755.330 |
|  12 | RID-GNSS-010  | 99.455426 | 14.491049 |   112 |      -8.0 | 45,194.856 | 52,682.966 |
|  13 | RID-GNSS-011  | 99.554406 | 14.312074 |   101 |      -6.2 | 55,869.703 | 32,880.767 |
|  14 | RID-GNSS-012  | 99.554651 | 14.310476 |    99 |      -5.9 | 55,896.207 | 32,703.976 |
|  15 | RID-GNSS-013  | 99.586083 | 14.380066 |    97 |      -4.9 | 59,284.484 | 40,404.598 |
|  16 | RID-GNSS-014  | 99.586888 | 14.378912 |    96 |      -4.7 | 59,371.298 | 40,276.965 |
|  17 | RID-GNSS-015  | 99.637282 | 14.428737 |    97 |      -3.2 | 64,803.287 | 45,792.480 |
|  18 | RID-GNSS-016  | 99.637328 | 14.425973 |    94 |      -2.7 | 64,808.401 | 45,486.582 |
|  19 | RID-GNSS-017  | 99.658843 | 14.498068 |    97 |      -2.3 | 67,122.908 | 53,465.001 |
|  20 | RID-GNSS-018  | 99.660703 | 14.497375 |    95 |      -1.9 | 67,323.485 | 53,388.453 |
|  21 | RID-GNSS-019  | 99.639670 | 14.588679 |    87 |      -1.5 | 65,050.026 | 63,489.347 |
|  22 | RID-GNSS-020  | 99.641275 | 14.587939 |    88 |      -1.6 | 65,222.957 | 63,407.523 |
|  23 | RID-GNSS-021  | 99.641666 | 14.644517 |    87 |      -1.5 | 65,261.235 | 69,667.725 |
|  24 | RID-GNSS-022  | 99.642754 | 14.643998 |    85 |      -1.1 | 65,378.467 | 69,610.369 |
|  25 | RID-GNSS-023  | 99.688102 | 14.727748 |    88 |       0.7 | 70,255.974 | 78,880.618 |
|  26 | RID-GNSS-024  | 99.689335 | 14.728310 |    86 |       0.9 | 70,388.660 | 78,942.845 |
|  27 | RID-GNSS-025  | 99.471661 | 14.484961 |   104 |      -6.9 | 46,944.962 | 52,009.092 |
|  28 | RID-GNSS-026  | 99.473050 | 14.484892 |   103 |      -6.8 | 47,094.648 | 52,001.427 |
|  29 | RID-GNSS-027  | 99.496875 | 14.489702 |   103 |      -6.8 | 49,663.069 | 52,533.495 |
|  30 | RID-GNSS-028  | 99.498339 | 14.490006 |   105 |      -7.2 | 49,820.950 | 52,567.093 |
|  31 | RID-GNSS-029  | 99.500047 | 14.458872 |    97 |      -5.9 | 50,005.069 | 49,122.302 |
|  32 | RID-GNSS-030  | 99.501365 | 14.459113 |    99 |      -6.2 | 50,147.138 | 49,148.950 |
|  33 | RID-GNSS-031  | 99.519507 | 14.328425 |    93 |      -5.2 | 52,104.415 | 34,689.332 |
|  34 | RID-GNSS-032  | 99.521162 | 14.329150 |    95 |      -5.6 | 52,282.945 | 34,769.507 |
|  35 | RID-GNSS-033  | 99.603353 | 14.401571 |   102 |      -5.1 | 61,146.047 | 42,784.843 |
|  36 | RID-GNSS-034  | 99.603098 | 14.400134 |   102 |      -5.2 | 61,118.644 | 42,625.821 |
|  37 | RID-GNSS-035  | 99.663234 | 14.526812 |    94 |      -1.7 | 67,593.993 | 56,645.740 |
|  38 | RID-GNSS-036  | 99.661689 | 14.527681 |    95 |      -1.9 | 67,427.422 | 56,741.733 |
|  39 | RID-GNSS-037  | 99.637304 | 14.567171 |    96 |      -3.1 | 64,796.480 | 61,109.454 |
|  40 | RID-GNSS-038  | 99.638886 | 14.567064 |    95 |      -2.9 | 64,966.940 | 61,097.714 |
|  41 | RID-GNSS-039  | 99.655834 | 14.654815 |    81 |       0.0 | 66,786.649 | 70,808.198 |
|  42 | RID-GNSS-040  | 99.655875 | 14.655976 |    82 |      -0.1 | 66,790.996 | 70,936.651 |
|  43 | RID-GNSS-041  | 99.678825 | 14.685007 |    87 |       0.2 | 69,260.720 | 74,150.589 |
|  44 | RID-GNSS-042  | 99.678732 | 14.686423 |    88 |       0.0 | 69,250.524 | 74,307.287 |
|  45 | RID-GNSS-043  | 99.688775 | 14.704124 |    88 |       0.6 | 70,330.613 | 76,266.770 |
|  46 | RID-GNSS-044  | 99.689977 | 14.703518 |    87 |       0.8 | 70,460.134 | 76,199.748 |
|  47 | RID-GNSS-045  | 99.166861 | 14.540621 |   203 |      -6.5 | 14,095.074 | 58,193.603 |
|  48 | RID-GNSS-046  | 99.165321 | 14.541526 |   195 |      -5.1 | 13,929.197 | 58,293.905 |
|  49 | RID-GNSS-047  | 99.364075 | 14.552361 |   170 |     -14.7 | 35,351.159 | 59,470.693 |
|  50 | RID-GNSS-048  | 99.363161 | 14.553795 |   172 |     -15.0 | 35,252.790 | 59,629.411 |
|  51 | RID-GNSS-049  | 99.687500 | 14.752835 |    83 |       1.3 | 70,188.848 | 81,656.347 |
|  52 | RID-GNSS-050  | 99.688426 | 14.754319 |    80 |       1.8 | 70,288.384 | 81,820.620 |
|  53 | RID-GNSS-051  | 99.728150 | 14.736188 |    53 |       8.4 | 74,567.632 | 79,818.470 |
|  54 | RID-GNSS-052  | 99.729709 | 14.735548 |    52 |       8.7 | 74,735.545 | 79,747.812 |
|  55 | RID-GNSS-053  | 99.707611 | 14.721819 |    69 |       4.6 | 72,357.383 | 78,226.439 |
|  56 | RID-GNSS-054  | 99.708233 | 14.724233 |    66 |       5.2 | 72,424.091 | 78,493.586 |
|  57 | RID-GNSS-055  | 99.733270 | 14.716708 |    52 |       9.0 | 75,121.216 | 77,663.572 |
|  58 | RID-GNSS-056  | 99.733287 | 14.718287 |    53 |       8.9 | 75,122.841 | 77,838.356 |
|  59 | RID-GNSS-057  | 99.705559 | 14.689230 |    71 |       4.2 | 72,139.745 | 74,620.377 |
|  60 | RID-GNSS-058  | 99.706318 | 14.691084 |    71 |       4.3 | 72,221.235 | 74,825.564 |
|  61 | RID-GNSS-059  | 99.744384 | 14.692168 |    47 |      10.4 | 76,320.984 | 74,949.608 |
|  62 | RID-GNSS-060  | 99.744219 | 14.690872 |    47 |      10.4 | 76,303.363 | 74,806.190 |
|  63 | RID-GNSS-061  | 99.720258 | 14.678788 |    61 |       6.7 | 73,723.968 | 73,466.466 |
|  64 | RID-GNSS-062  | 99.721534 | 14.678697 |    60 |       6.9 | 73,861.444 | 73,456.501 |
|  65 | RID-GNSS-063  | 99.757100 | 14.678713 |    44 |      11.9 | 77,692.283 | 73,462.383 |
|  66 | RID-GNSS-064  | 99.756907 | 14.677442 |    44 |      11.8 | 77,671.701 | 73,321.690 |
|  67 | RID-GNSS-065  | 99.716027 | 14.662183 |    60 |       6.5 | 73,270.029 | 71,628.695 |
|  68 | RID-GNSS-066  | 99.716230 | 14.660308 |    59 |       6.7 | 73,292.085 | 71,421.266 |
|  69 | RID-GNSS-067  | 99.762514 | 14.654951 |    42 |      12.5 | 78,278.426 | 70,833.810 |
|  70 | RID-GNSS-068  | 99.761800 | 14.652950 |    43 |      12.4 | 78,201.831 | 70,612.369 |
|  71 | RID-GNSS-069  | 99.815097 | 14.640681 |    31 |      18.7 | 83,945.058 | 69,262.106 |
|  72 | RID-GNSS-070  | 99.814664 | 14.638525 |    31 |      18.7 | 83,898.748 | 69,023.508 |
|  73 | RID-GNSS-071  | 99.694967 | 14.649930 |    64 |       4.7 | 71,002.600 | 70,270.873 |
|  74 | RID-GNSS-072  | 99.694978 | 14.648048 |    64 |       4.7 | 71,004.042 | 70,062.642 |
|  75 | RID-GNSS-073  | 99.745258 | 14.627599 |    45 |      10.8 | 76,422.893 | 67,805.357 |
|  76 | RID-GNSS-074  | 99.743296 | 14.627496 |    46 |      10.5 | 76,211.504 | 67,793.690 |
|  77 | RID-GNSS-075  | 99.796041 | 14.632125 |    33 |      16.7 | 81,893.346 | 68,312.592 |
|  78 | RID-GNSS-076  | 99.795599 | 14.630569 |    33 |      16.7 | 81,845.983 | 68,140.390 |
|  79 | RID-GNSS-077  | 99.795224 | 14.604054 |    32 |      16.8 | 81,809.410 | 65,206.558 |
|  80 | RID-GNSS-078  | 99.796561 | 14.604689 |    32 |      16.9 | 81,953.374 | 65,276.954 |
|  81 | RID-GNSS-079  | 99.686829 | 14.623462 |    71 |       3.1 | 70,128.354 | 67,341.590 |
|  82 | RID-GNSS-080  | 99.686176 | 14.622080 |    72 |       3.0 | 70,058.176 | 67,188.592 |
|  83 | RID-GNSS-081  | 99.733082 | 14.611961 |    45 |      10.0 | 75,112.827 | 66,073.693 |
|  84 | RID-GNSS-082  | 99.732725 | 14.610588 |    45 |      10.0 | 75,074.587 | 65,921.725 |
|  85 | RID-GNSS-083  | 99.687746 | 14.605499 |    57 |       5.4 | 70,228.818 | 65,354.109 |
|  86 | RID-GNSS-084  | 99.687165 | 14.604063 |    57 |       5.3 | 70,166.370 | 65,195.238 |
|  87 | RID-GNSS-085  | 99.668020 | 14.581001 |    78 |       1.1 | 68,105.426 | 62,641.841 |
|  88 | RID-GNSS-086  | 99.669312 | 14.580351 |    77 |       1.4 | 68,244.744 | 62,570.013 |
|  89 | RID-GNSS-087  | 99.693208 | 14.586525 |    56 |       5.9 | 70,819.159 | 63,255.233 |
|  90 | RID-GNSS-088  | 99.694168 | 14.585421 |    56 |       6.0 | 70,922.655 | 63,133.116 |
|  91 | RID-GNSS-089  | 99.731174 | 14.600769 |    45 |       9.9 | 74,908.553 | 64,835.056 |
|  92 | RID-GNSS-090  | 99.732201 | 14.600850 |    44 |      10.1 | 75,019.239 | 64,844.120 |
|  93 | RID-GNSS-091  | 99.670141 | 14.574337 |    68 |       2.7 | 68,334.537 | 61,904.714 |
|  94 | RID-GNSS-092  | 99.671937 | 14.573794 |    67 |       3.1 | 68,528.147 | 61,844.720 |
|  95 | RID-GNSS-093  | 99.658608 | 14.553350 |    85 |      -0.4 | 67,093.329 | 59,581.687 |
|  96 | RID-GNSS-094  | 99.660682 | 14.554499 |    82 |       0.1 | 67,316.755 | 59,709.000 |
|  97 | RID-GNSS-095  | 99.710225 | 14.573544 |    59 |       6.4 | 72,654.166 | 61,820.517 |
|  98 | RID-GNSS-096  | 99.713118 | 14.571322 |    66 |       5.5 | 72,966.096 | 61,574.997 |
|  99 | RID-GNSS-097  | 99.749058 | 14.600453 |    41 |      11.7 | 76,835.544 | 64,802.172 |
| 100 | RID-GNSS-098  | 99.751505 | 14.599493 |    40 |      12.0 | 77,099.359 | 64,696.232 |
| 101 | RID-GNSS-099  | 99.701203 | 14.543222 |    67 |       4.5 | 71,684.835 | 58,464.701 |
| 102 | RID-GNSS-100  | 99.702252 | 14.542080 |    68 |       4.6 | 71,798.013 | 58,338.388 |
| 103 | RID-GNSS-101  | 99.737489 | 14.576657 |    56 |       8.6 | 75,591.749 | 62,167.944 |
| 104 | RID-GNSS-102  | 99.737345 | 14.577737 |    55 |       8.7 | 75,576.181 | 62,287.407 |
| 105 | RID-GNSS-103  | 99.760461 | 14.597793 |    40 |      12.8 | 78,064.540 | 64,509.178 |
| 106 | RID-GNSS-104  | 99.761524 | 14.598111 |    41 |      12.7 | 78,179.044 | 64,544.539 |
| 107 | RID-GNSS-105  | 99.766607 | 14.578337 |    40 |      13.3 | 78,729.345 | 62,357.219 |
| 108 | RID-GNSS-106  | 99.767685 | 14.579784 |    40 |      13.2 | 78,845.360 | 62,517.493 |
| 109 | RID-GNSS-107  | 99.768618 | 14.548066 |    50 |      11.7 | 78,949.978 | 59,008.115 |
| 110 | RID-GNSS-108  | 99.768535 | 14.549066 |    49 |      11.9 | 78,940.902 | 59,118.728 |
| 111 | RID-GNSS-109  | 99.723825 | 14.525301 |    64 |       6.5 | 74,124.944 | 56,484.093 |
| 112 | RID-GNSS-110  | 99.725094 | 14.523757 |    65 |       6.4 | 74,261.840 | 56,313.383 |
| 113 | RID-GNSS-111  | 99.762054 | 14.523128 |    60 |       9.7 | 78,245.722 | 56,247.996 |
| 114 | RID-GNSS-112  | 99.762451 | 14.525825 |    59 |       9.9 | 78,288.225 | 56,546.453 |
| 115 | RID-GNSS-113  | 99.761134 | 14.509647 |    67 |       8.5 | 78,148.259 | 54,756.293 |
| 116 | RID-GNSS-114  | 99.763471 | 14.508717 |    67 |       8.7 | 78,400.281 | 54,653.751 |
| 117 | RID-GNSS-115  | 99.789258 | 14.505498 |    73 |       9.8 | 81,180.464 | 54,300.871 |
| 118 | RID-GNSS-116  | 99.789027 | 14.503554 |    74 |       9.6 | 81,155.781 | 54,085.783 |
| 119 | RID-GNSS-117  | 99.709484 | 14.514695 |    68 |       5.0 | 72,580.242 | 55,309.140 |
| 120 | RID-GNSS-118  | 99.709453 | 14.516464 |    68 |       5.0 | 72,576.781 | 55,504.900 |
| 121 | RID-GNSS-119  | 99.691728 | 14.498472 |    74 |       3.0 | 70,667.797 | 53,512.423 |
| 122 | RID-GNSS-120  | 99.690255 | 14.499156 |    76 |       2.6 | 70,508.965 | 53,587.965 |
| 123 | RID-GNSS-121  | 99.723370 | 14.492606 |    74 |       4.8 | 74,079.393 | 52,866.534 |
| 124 | RID-GNSS-122  | 99.725132 | 14.492070 |    72 |       5.3 | 74,269.373 | 52,807.386 |
| 125 | RID-GNSS-123  | 99.774028 | 14.489141 |    71 |       9.0 | 79,540.878 | 52,489.011 |
| 126 | RID-GNSS-124  | 99.774295 | 14.487756 |    71 |       8.9 | 79,569.820 | 52,335.860 |
| 127 | RID-GNSS-125  | 99.693078 | 14.477380 |    69 |       3.8 | 70,815.358 | 51,178.872 |
| 128 | RID-GNSS-126  | 99.692867 | 14.475732 |    68 |       3.9 | 70,792.784 | 50,996.478 |
| 129 | RID-GNSS-127  | 99.714963 | 14.426613 |    42 |       9.4 | 73,180.002 | 45,563.872 |
| 130 | RID-GNSS-128  | 99.716913 | 14.426584 |    41 |       9.6 | 73,390.239 | 45,560.848 |
| 131 | RID-GNSS-129  | 99.750341 | 14.466167 |    57 |       9.4 | 76,990.071 | 49,944.185 |
| 132 | RID-GNSS-130  | 99.749555 | 14.464015 |    56 |       9.4 | 76,905.622 | 49,705.985 |
| 133 | RID-GNSS-131  | 99.752495 | 14.408349 |    39 |      12.3 | 77,229.401 | 43,547.123 |
| 134 | RID-GNSS-132  | 99.752163 | 14.406261 |    39 |      12.4 | 77,193.842 | 43,316.142 |
| 135 | RID-GNSS-133  | 99.674629 | 14.436152 |    68 |       3.0 | 68,829.808 | 46,615.594 |
| 136 | RID-GNSS-134  | 99.674669 | 14.434049 |    72 |       2.3 | 68,834.333 | 46,382.884 |
| 137 | RID-GNSS-135  | 99.719665 | 14.405527 |    35 |      10.6 | 73,689.244 | 43,231.312 |
| 138 | RID-GNSS-136  | 99.718165 | 14.406385 |    36 |      10.4 | 73,527.384 | 43,326.111 |
| 139 | RID-GNSS-137  | 99.667409 | 14.409900 |    61 |       3.8 | 68,053.482 | 43,710.353 |
| 140 | RID-GNSS-138  | 99.669360 | 14.409892 |    64 |       3.4 | 68,263.863 | 43,709.689 |
| 141 | RID-GNSS-139  | 99.684880 | 14.380455 |    41 |       7.7 | 69,940.167 | 40,453.883 |
| 142 | RID-GNSS-140  | 99.686511 | 14.379477 |    40 |       7.9 | 70,116.162 | 40,345.872 |
| 143 | RID-GNSS-141  | 99.652371 | 14.415809 |    72 |       1.4 | 66,431.299 | 44,363.079 |
| 144 | RID-GNSS-142  | 99.652317 | 14.414112 |    70 |       1.7 | 66,425.555 | 44,175.316 |
| 145 | RID-GNSS-143  | 99.669022 | 14.381103 |    49 |       5.8 | 68,229.713 | 40,524.372 |
| 146 | RID-GNSS-144  | 99.670737 | 14.380113 |    47 |       6.0 | 68,414.756 | 40,414.874 |
| 147 | RID-GNSS-145  | 99.644939 | 14.397687 |    69 |       1.4 | 65,631.068 | 42,357.475 |
| 148 | RID-GNSS-146  | 99.645064 | 14.396117 |    68 |       1.7 | 65,644.751 | 42,183.747 |
| 149 | RID-GNSS-147  | 99.681027 | 14.360542 |    41 |       7.6 | 69,526.337 | 38,250.337 |
| 150 | RID-GNSS-148  | 99.682250 | 14.360305 |    42 |       7.5 | 69,658.257 | 38,224.289 |
| 151 | RID-GNSS-149  | 99.621412 | 14.373073 |    75 |      -0.4 | 63,095.208 | 39,632.611 |
| 152 | RID-GNSS-150  | 99.621746 | 14.372025 |    77 |      -0.7 | 63,131.333 | 39,516.679 |
| 153 | RID-GNSS-151  | 99.638069 | 14.358826 |    70 |       1.0 | 64,892.757 | 38,057.300 |
| 154 | RID-GNSS-152  | 99.638122 | 14.357459 |    69 |       1.1 | 64,898.627 | 37,906.016 |
| 155 | RID-GNSS-153  | 99.595585 | 14.362570 |    80 |      -2.0 | 60,310.053 | 38,469.232 |
| 156 | RID-GNSS-154  | 99.594318 | 14.362464 |    81 |      -2.1 | 60,173.406 | 38,457.409 |
| 157 | RID-GNSS-155  | 99.635049 | 14.339189 |    56 |       3.1 | 64,568.290 | 35,884.362 |
| 158 | RID-GNSS-156  | 99.637192 | 14.338243 |    57 |       3.0 | 64,799.579 | 35,779.853 |
| 159 | RID-GNSS-157  | 99.609680 | 14.326041 |    66 |       0.6 | 61,832.302 | 34,428.210 |
| 160 | RID-GNSS-158  | 99.611089 | 14.325989 |    65 |       0.8 | 61,984.282 | 34,422.503 |
| 161 | RID-GNSS-159  | 99.582010 | 14.337042 |    88 |      -3.6 | 58,846.872 | 35,644.136 |
| 162 | RID-GNSS-160  | 99.583411 | 14.338557 |    86 |      -3.2 | 58,997.893 | 35,811.823 |
| 163 | RID-GNSS-161  | 99.600523 | 14.311107 |    85 |      -2.6 | 60,845.147 | 32,775.478 |
| 164 | RID-GNSS-162  | 99.598708 | 14.311169 |    88 |      -3.1 | 60,649.328 | 32,782.201 |
| 165 | RID-GNSS-163  | 99.665222 | 14.299128 |    39 |       7.0 | 67,826.335 | 31,454.022 |
| 166 | RID-GNSS-164  | 99.664828 | 14.297822 |    39 |       7.1 | 67,783.985 | 31,309.543 |
| 167 | RID-GNSS-165  | 99.608313 | 14.294164 |    76 |      -1.0 | 61,686.456 | 30,901.182 |
| 168 | RID-GNSS-166  | 99.608019 | 14.293116 |    79 |      -1.4 | 61,654.803 | 30,785.242 |
| 169 | RID-GNSS-167  | 99.635288 | 14.294844 |    51 |       3.9 | 64,596.985 | 30,977.975 |
| 170 | RID-GNSS-168  | 99.636712 | 14.294827 |    50 |       4.1 | 64,750.624 | 30,976.130 |
| 171 | RID-GNSS-169  | 99.539682 | 14.290166 |    95 |      -5.5 | 54,281.545 | 30,456.499 |
| 172 | RID-GNSS-170  | 99.539020 | 14.288440 |    95 |      -5.3 | 54,210.172 | 30,265.508 |
| 173 | RID-GNSS-171  | 99.543372 | 14.254574 |    75 |      -2.2 | 54,680.413 | 26,518.535 |
| 174 | RID-GNSS-172  | 99.544038 | 14.252100 |    77 |      -2.6 | 54,752.373 | 26,244.807 |
| 175 | RID-GNSS-173  | 99.521159 | 14.263374 |    78 |      -2.8 | 52,283.284 | 27,491.898 |
| 176 | RID-GNSS-174  | 99.522098 | 14.261514 |    75 |      -2.5 | 52,384.635 | 27,286.154 |
| 177 | RID-GNSS-175  | 99.520699 | 14.213118 |    71 |      -1.8 | 52,234.115 | 21,931.467 |
| 178 | RID-GNSS-176  | 99.521716 | 14.213376 |    72 |      -1.9 | 52,343.896 | 21,960.045 |



