# RTK_Construct
```sh
usage: constr_LDP [-h] [-c] [-u] [-o OFFSET_PP] LDP_TOML

Calculate CSF from points defined by GPKG, CSV or single point

positional arguments:
  LDP_TOML

options:
  -h, --help            show this help message and exit
  -c, --csf             show CSF table
  -u, --utm             show UTM-LDP table
  -o OFFSET_PP, --OFFSET_PP OFFSET_PP
                        offset for project plane
  
  ตัวอย่างการออกแบบระบบพิกัดความคลาดเคลื่อนต่ำสำหรับพื้นที่โครงการ “บ่อพลอย” ได้จากการพิจารณากำหนดจุดสุ่มทดสอบใช้หมุดควบคุม RID-GNSS-001 ถึง RID-GNSS-176 ทั้งหมดในโครงการ พร้อมจุดตัวแทนแนวกึ่งกลางของอุโมงค์ส่งน้ำ T1, T2 และ T3 กำหนดค่าออกแบบเบื้องต้นเป็น  

|     | Point         |     UTM_N |    UTM_E  |    MSL  |
|----:|:--------------|----------:|----------:|--------:|
|   1 | Tunnel-Inlet  | 1,608,603 |  518,123  |   160.0 |
|   2 | Tunnel-Mid    | 1,608,118 |  529,367  |   158.5 |
|   3 | Tunnel-Outlet | 1,607,439 |  539,639  |   157.0 |

  ความสูงเหนือทะเลปานกลางของหมุดควบคุมได้จากงานเดินระดับดิฟเฟอเรนเชียลงานชั้น 2 โครงการเลือกใช้การฉายแผนที่ระบบพิกัดทรานเวอร์สเมอร์เคเตอร์ (Transverse Mercator :TM ) ให้ศูนย์กลางเมอริเดียนเป็น 99°30’ การคำนวณไปความสูงเหนือทรงรีใช้แบบจำลองยีออยด์ TGM-2017 เนื่องจากมีงานออกแบบและก่อสร้างอยูลึกใต้ผิวภูมิประเทศลงไป จึงมีการกำหนดจุดทดสอบ T1, T2 และ T3 ลึกใต้ผิวภูมิประเทศ ในขั้นต้นระดับ projection plane อยู่ในเกณฑ์ค่าเฉลี่ยของระดับทะเลปานกลางของจุดทดสอบ (offset PP=0 m) ได้ค่าเสกลแฟกเตอร์ร่วม Points CSF min/mean/max [ppm] : -7.4 / +2.8 / +18.4 ซึ่งไม่เกิน +20 ppm ดังนั้นเพื่อให้ระนาบการฉายแผนที่อยู่ในเกณฑ์ค่าระดับของอุโมงค์ด้วยจึงเปลี่ยนระดับการฉายต่ำลง offset PP = -30 m และ CSF ลดต่ำลงไปอีกทั้งในทิศทางเหนือระนาบและใต้ระนาบลึกลงไป ดั้งนั้นตัดสินใจเลือกระนาบและได้ ค่าเสกลแฟกเตอร์ร่วม Points CSF min/mean/max [ppm] : -12.4 / -2.2 / +13.4

## แผนที่แสดงที่ตั้งจุดตัวแทนของโครงการก่อสร้างเพื่อใช้คำนวนสเกลแฟกเตอร์ร่วม (Combined Scale Factor : CSF )
![Map](https://github.com/phisan-chula/RTK_Construct/blob/main/BoPloi_LDP.png)

## ช่วงของค่าสเกลแฟกเตอร์ร่วม (Combined Scale Factor : CSF ) ได้จากการกำหนดระนาบการฉาย TM ทีระนาบการฉาย ณ ระนาบก่อสร้าง ณ รทก ต่างๆ  

Offset Projection Plane = +20.0 m. <defined by CLI args>  
Designed Project Plane : hPP=+62.54 m. / MSL=+96.39 m. tied to k0 = 1.000010  
Points CSF min/mean/max [ppm] : -4.4 / +5.8 / +21.4  

Offset Projection Plane = 0.0 m. <defined by CLI args>  
Designed Project Plane : hPP=+42.54 m. / MSL=+76.39 m. tied to k0 = 1.000007  
Points CSF min/mean/max [ppm] : -7.4 / +2.8 / +18.4  

Offset Projection Plane = -20.0 m. <defined by CLI args>  
Designed Project Plane : hPP=+22.54 m. / MSL=+56.39 m. tied to k0 = 1.000004  
Points CSF min/mean/max [ppm] : -10.4 / -0.2 / +15.4  

Offset Projection Plane = -30.0 m. <defined by CLI args>  
Designed Project Plane : hPP=+12.54 m. / MSL=+46.39 m. tied to k0 = 1.000002  
Points CSF min/mean/max [ppm] : -12.4 / -2.2 / +13.4  

Offset Projection Plane = -40.0 m. <defined by CLI args>  
Designed Project Plane : hPP=+2.54 m. / MSL=+36.39 m. tied to k0 = 1.000000  
Points CSF min/mean/max [ppm] : -14.4 / -4.2 / +11.4   

ทำการปรับค่า False Easting และ Fase Northing ให้พิกัดท้องถิ่นมีตัวเลขกระทัดรัดและกระชับ ไม่มีตัวเลขค่าพิกัดมีเครื่องหมายลบ ในมาตรฐานฉบับนี้จึงเสนอว่า ค่าพารามิเตอร์ระบบพิกัดท้องถิ่น LDP ของโครงการ “บ่อพลอย” (LDP-BOPLOI) คือ 

## ค่าพารามิเตอร์ระบบพิกัดท้องถิ่น LDP ของโครงการ “บ่อพลอย” (LDP-BOPLOI)  
  
================================ LDP Defintion ===============================  
+proj=tmerc +lat_0=0.0 +lon_0=99.63333333333334 +k_0=1.000002  +x_0=60000  +y_0=-1550000  
        +a=6378137.0 +b=6356752.3142 +units=m +no_defs +type=crs  
+proj=tmerc +lat_0=0.0 +lon_0=99°38′ +k_0=1.000002  +x_0=60000  +y_0=-1550000  
        +a=6378137.0 +b=6356752.3142 +units=m +no_defs +type=crs  


## ตาราง แสดงการคำนวน LDP และค่า LDP ที่เกิดขึ้นเหมาะการนำไปใช้งาน สดวกและแม่นยำ

|     | Point         |       lng |       lat |    MSL |   CSF_ppm |      LDP_E |      LDP_N |
|----:|:--------------|----------:|----------:|-------:|----------:|-----------:|-----------:|
|   0 | Tunnel-Inlet  | 99.168228 | 14.550374 | +160.0 |     +13.4 |  9,874.091 | 59,294.344 |
|   1 | Tunnel-Mid    | 99.272595 | 14.545890 | +158.5 |      +1.1 | 21,121.419 | 58,777.909 |
|   2 | Tunnel-Outlet | 99.367932 | 14.539621 | +157.0 |      -7.2 | 31,395.619 | 58,070.116 |
|   3 | RID-GNSS-001  | 99.168657 | 14.552301 | +185.3 |      +9.3 |  9,920.767 | 59,507.498 |
|   4 | RID-GNSS-002  | 99.168475 | 14.551096 | +187.5 |      +9.0 |  9,900.939 | 59,374.187 |
|   5 | RID-GNSS-003  | 99.360889 | 14.546048 | +192.8 |     -12.3 | 30,637.369 | 58,782.166 |
|   6 | RID-GNSS-004  | 99.359894 | 14.545350 | +194.0 |     -12.4 | 30,530.102 | 58,704.979 |
|   7 | RID-GNSS-005  | 99.375530 | 14.539136 | +160.6 |      -8.3 | 32,214.540 | 58,015.537 |
|   8 | RID-GNSS-006  | 99.374453 | 14.537768 | +157.5 |      -7.8 | 32,098.278 | 57,864.237 |
|   9 | RID-GNSS-007  | 99.420384 | 14.521585 | +133.9 |      -7.2 | 37,046.968 | 56,068.583 |
|  10 | RID-GNSS-008  | 99.420039 | 14.520434 | +131.8 |      -6.8 | 37,009.694 | 55,941.328 |
|  11 | RID-GNSS-009  | 99.454297 | 14.491703 | +113.7 |      -5.9 | 40,699.808 | 52,759.179 |
|  12 | RID-GNSS-010  | 99.455426 | 14.491049 | +112.4 |      -5.7 | 40,821.390 | 52,686.744 |
|  13 | RID-GNSS-011  | 99.554406 | 14.312074 | +101.2 |      -7.7 | 51,484.778 | 32,878.362 |
|  14 | RID-GNSS-012  | 99.554651 | 14.310476 |  +99.4 |      -7.4 | 51,511.180 | 32,701.556 |
|  15 | RID-GNSS-013  | 99.586083 | 14.380066 |  +97.4 |      -7.7 | 54,903.889 | 40,400.207 |
|  16 | RID-GNSS-014  | 99.586888 | 14.378912 |  +96.1 |      -7.5 | 54,990.629 | 40,272.525 |
|  17 | RID-GNSS-015  | 99.637282 | 14.428737 |  +96.7 |      -7.9 | 60,425.789 | 45,784.873 |
|  18 | RID-GNSS-016  | 99.637328 | 14.425973 |  +93.5 |      -7.4 | 60,430.727 | 45,478.974 |
|  19 | RID-GNSS-017  | 99.658843 | 14.498068 |  +96.8 |      -7.8 | 62,749.858 | 53,456.006 |
|  20 | RID-GNSS-018  | 99.660703 | 14.497375 |  +95.1 |      -7.5 | 62,950.389 | 53,379.341 |
|  21 | RID-GNSS-019  | 99.639670 | 14.588679 |  +87.0 |      -6.3 | 60,682.845 | 63,481.509 |
|  22 | RID-GNSS-020  | 99.641275 | 14.587939 |  +87.9 |      -6.5 | 60,855.727 | 63,399.584 |
|  23 | RID-GNSS-021  | 99.641666 | 14.644517 |  +87.2 |      -6.4 | 60,897.681 | 69,659.732 |
|  24 | RID-GNSS-022  | 99.642754 | 14.643998 |  +85.2 |      -6.0 | 61,014.878 | 69,602.308 |
|  25 | RID-GNSS-023  | 99.688102 | 14.727748 |  +87.6 |      -6.0 | 65,897.826 | 78,869.624 |
|  26 | RID-GNSS-024  | 99.689335 | 14.728310 |  +86.2 |      -5.8 | 66,030.547 | 78,931.772 |
|  27 | RID-GNSS-025  | 99.471661 | 14.484961 | +104.2 |      -5.3 | 42,571.107 | 52,011.849 |
|  28 | RID-GNSS-026  | 99.473050 | 14.484892 | +103.3 |      -5.2 | 42,720.789 | 52,004.098 |
|  29 | RID-GNSS-027  | 99.496875 | 14.489702 | +102.9 |      -6.1 | 45,289.522 | 52,534.671 |
|  30 | RID-GNSS-028  | 99.498339 | 14.490006 | +105.4 |      -6.6 | 45,447.423 | 52,568.177 |
|  31 | RID-GNSS-029  | 99.500047 | 14.458872 |  +97.1 |      -5.4 | 45,629.538 | 49,123.277 |
|  32 | RID-GNSS-030  | 99.501365 | 14.459113 |  +98.8 |      -5.7 | 45,771.622 | 49,149.843 |
|  33 | RID-GNSS-031  | 99.519507 | 14.328425 |  +92.6 |      -5.4 | 47,720.535 | 34,689.092 |
|  34 | RID-GNSS-032  | 99.521162 | 14.329150 |  +95.4 |      -5.9 | 47,899.111 | 34,769.164 |
|  35 | RID-GNSS-033  | 99.603353 | 14.401571 | +101.8 |      -8.6 | 56,766.823 | 42,779.368 |
|  36 | RID-GNSS-034  | 99.603098 | 14.400134 | +102.1 |      -8.6 | 56,739.328 | 42,620.362 |
|  37 | RID-GNSS-035  | 99.663234 | 14.526812 |  +94.4 |      -7.4 | 63,222.795 | 56,636.452 |
|  38 | RID-GNSS-036  | 99.661689 | 14.527681 |  +95.4 |      -7.6 | 63,056.281 | 56,732.541 |
|  39 | RID-GNSS-037  | 99.637304 | 14.567171 |  +96.3 |      -7.8 | 60,427.906 | 61,101.776 |
|  40 | RID-GNSS-038  | 99.638886 | 14.567064 |  +95.3 |      -7.6 | 60,598.358 | 61,089.937 |
|  41 | RID-GNSS-039  | 99.655834 | 14.654815 |  +81.4 |      -5.4 | 62,423.758 | 70,799.302 |
|  42 | RID-GNSS-040  | 99.655875 | 14.655976 |  +82.4 |      -5.5 | 62,428.181 | 70,927.752 |
|  43 | RID-GNSS-041  | 99.678825 | 14.685007 |  +87.1 |      -6.1 | 64,899.784 | 74,140.214 |
|  44 | RID-GNSS-042  | 99.678732 | 14.686423 |  +88.3 |      -6.2 | 64,889.681 | 74,296.917 |
|  45 | RID-GNSS-043  | 99.688775 | 14.704124 |  +88.0 |      -6.1 | 65,970.919 | 76,255.750 |
|  46 | RID-GNSS-044  | 99.689977 | 14.703518 |  +87.1 |      -5.9 | 66,100.400 | 76,188.651 |
|  47 | RID-GNSS-045  | 99.166861 | 14.540621 | +202.6 |      +6.8 |  9,724.586 | 58,215.563 |
|  48 | RID-GNSS-046  | 99.165321 | 14.541526 | +194.9 |      +8.3 |  9,558.765 | 58,315.964 |
|  49 | RID-GNSS-047  | 99.364075 | 14.552361 | +170.4 |      -9.0 | 30,981.616 | 59,480.241 |
|  50 | RID-GNSS-048  | 99.363161 | 14.553795 | +172.1 |      -9.2 | 30,883.339 | 59,639.018 |
|  51 | RID-GNSS-049  | 99.687500 | 14.752835 |  +83.3 |      -5.3 | 65,832.344 | 81,645.374 |
|  52 | RID-GNSS-050  | 99.688426 | 14.754319 |  +80.5 |      -4.9 | 65,931.976 | 81,809.587 |
|  53 | RID-GNSS-051  | 99.728150 | 14.736188 |  +53.4 |      +0.2 | 70,210.006 | 79,804.917 |
|  54 | RID-GNSS-052  | 99.729709 | 14.735548 |  +52.5 |      +0.4 | 70,377.876 | 79,734.161 |
|  55 | RID-GNSS-053  | 99.707611 | 14.721819 |  +69.4 |      -2.8 | 67,998.832 | 78,214.207 |
|  56 | RID-GNSS-054  | 99.708233 | 14.724233 |  +66.0 |      -2.2 | 68,065.698 | 78,481.312 |
|  57 | RID-GNSS-055  | 99.733270 | 14.716708 |  +51.9 |      +0.6 | 70,762.311 | 77,649.710 |
|  58 | RID-GNSS-056  | 99.733287 | 14.718287 |  +52.5 |      +0.5 | 70,764.039 | 77,824.492 |
|  59 | RID-GNSS-057  | 99.705559 | 14.689230 |  +71.3 |      -3.1 | 67,779.066 | 74,608.300 |
|  60 | RID-GNSS-058  | 99.706318 | 14.691084 |  +70.9 |      -3.0 | 67,860.677 | 74,813.438 |
|  61 | RID-GNSS-059  | 99.744384 | 14.692168 |  +47.3 |      +1.6 | 71,960.465 | 74,935.061 |
|  62 | RID-GNSS-060  | 99.744219 | 14.690872 |  +47.3 |      +1.7 | 71,942.760 | 74,791.654 |
|  63 | RID-GNSS-061  | 99.720258 | 14.678788 |  +60.7 |      -1.1 | 69,362.596 | 73,453.463 |
|  64 | RID-GNSS-062  | 99.721534 | 14.678697 |  +60.0 |      -1.0 | 69,500.065 | 73,443.418 |
|  65 | RID-GNSS-063  | 99.757100 | 14.678713 |  +43.7 |      +2.6 | 73,330.874 | 73,447.041 |
|  66 | RID-GNSS-064  | 99.756907 | 14.677442 |  +44.1 |      +2.6 | 73,310.209 | 73,306.361 |
|  67 | RID-GNSS-065  | 99.716027 | 14.662183 |  +60.4 |      -1.2 | 68,907.578 | 71,615.975 |
|  68 | RID-GNSS-066  | 99.716230 | 14.660308 |  +59.4 |      -1.0 | 68,929.511 | 71,408.534 |
|  69 | RID-GNSS-067  | 99.762514 | 14.654951 |  +42.4 |      +3.0 | 73,915.462 | 70,818.148 |
|  70 | RID-GNSS-068  | 99.761800 | 14.652950 |  +43.0 |      +2.9 | 73,838.738 | 70,596.753 |
|  71 | RID-GNSS-069  | 99.815097 | 14.640681 |  +30.8 |      +7.2 | 79,581.110 | 69,243.125 |
|  72 | RID-GNSS-070  | 99.814664 | 14.638525 |  +30.7 |      +7.2 | 79,534.660 | 69,004.558 |
|  73 | RID-GNSS-071  | 99.694967 | 14.649930 |  +64.2 |      -2.2 | 66,639.366 | 70,259.498 |
|  74 | RID-GNSS-072  | 99.694978 | 14.648048 |  +64.2 |      -2.2 | 66,640.685 | 70,051.268 |
|  75 | RID-GNSS-073  | 99.745258 | 14.627599 |  +45.4 |      +2.0 | 72,058.166 | 67,790.814 |
|  76 | RID-GNSS-074  | 99.743296 | 14.627496 |  +46.2 |      +1.8 | 71,846.771 | 67,779.271 |
|  77 | RID-GNSS-075  | 99.796041 | 14.632125 |  +33.0 |      +5.9 | 77,528.862 | 68,294.829 |
|  78 | RID-GNSS-076  | 99.795599 | 14.630569 |  +32.8 |      +5.9 | 77,481.398 | 68,122.656 |
|  79 | RID-GNSS-077  | 99.795224 | 14.604054 |  +32.0 |      +6.0 | 77,443.103 | 65,188.878 |
|  80 | RID-GNSS-078  | 99.796561 | 14.604689 |  +31.6 |      +6.1 | 77,587.107 | 65,259.189 |
|  81 | RID-GNSS-079  | 99.686829 | 14.623462 |  +71.5 |      -3.5 | 65,763.404 | 67,330.749 |
|  82 | RID-GNSS-080  | 99.686176 | 14.622080 |  +72.2 |      -3.6 | 65,693.136 | 67,177.794 |
|  83 | RID-GNSS-081  | 99.733082 | 14.611961 |  +45.3 |      +1.6 | 70,747.094 | 66,059.935 |
|  84 | RID-GNSS-082  | 99.732725 | 14.610588 |  +45.0 |      +1.7 | 70,708.765 | 65,907.990 |
|  85 | RID-GNSS-083  | 99.687746 | 14.605499 |  +57.2 |      -1.2 | 65,862.700 | 65,343.223 |
|  86 | RID-GNSS-084  | 99.687165 | 14.604063 |  +57.4 |      -1.3 | 65,800.159 | 65,184.390 |
|  87 | RID-GNSS-085  | 99.668020 | 14.581001 |  +78.3 |      -4.8 | 63,737.732 | 62,632.218 |
|  88 | RID-GNSS-086  | 99.669312 | 14.580351 |  +76.8 |      -4.6 | 63,877.006 | 62,560.308 |
|  89 | RID-GNSS-087  | 99.693208 | 14.586525 |  +56.0 |      -1.0 | 66,451.806 | 63,244.016 |
|  90 | RID-GNSS-088  | 99.694168 | 14.585421 |  +55.8 |      -0.9 | 66,555.230 | 63,121.839 |
|  91 | RID-GNSS-089  | 99.731174 | 14.600769 |  +44.9 |      +1.6 | 70,542.095 | 64,821.428 |
|  92 | RID-GNSS-090  | 99.732201 | 14.600850 |  +44.1 |      +1.8 | 70,652.785 | 64,830.427 |
|  93 | RID-GNSS-091  | 99.670141 | 14.574337 |  +68.5 |      -3.2 | 63,966.409 | 61,894.961 |
|  94 | RID-GNSS-092  | 99.671937 | 14.573794 |  +66.6 |      -2.9 | 64,159.983 | 61,834.854 |
|  95 | RID-GNSS-093  | 99.658608 | 14.553350 |  +84.6 |      -5.9 | 62,723.850 | 59,572.674 |
|  96 | RID-GNSS-094  | 99.660682 | 14.554499 |  +82.0 |      -5.5 | 62,947.348 | 59,699.856 |
|  97 | RID-GNSS-095  | 99.710225 | 14.573544 |  +58.8 |      -1.1 | 68,285.959 | 61,808.235 |
|  98 | RID-GNSS-096  | 99.713118 | 14.571322 |  +66.0 |      -2.1 | 68,597.743 | 61,562.534 |
|  99 | RID-GNSS-097  | 99.749058 | 14.600453 |  +41.1 |      +2.8 | 72,469.049 | 64,787.414 |
| 100 | RID-GNSS-098  | 99.751505 | 14.599493 |  +40.5 |      +2.9 | 72,732.800 | 64,681.320 |
| 101 | RID-GNSS-099  | 99.701203 | 14.543222 |  +67.4 |      -2.6 | 67,314.673 | 58,453.011 |
| 102 | RID-GNSS-100  | 99.702252 | 14.542080 |  +67.5 |      -2.6 | 67,427.776 | 58,326.633 |
| 103 | RID-GNSS-101  | 99.737489 | 14.576657 |  +55.9 |      +0.1 | 71,223.721 | 62,153.938 |
| 104 | RID-GNSS-102  | 99.737345 | 14.577737 |  +55.1 |      +0.2 | 71,208.223 | 62,273.410 |
| 105 | RID-GNSS-103  | 99.760461 | 14.597793 |  +39.6 |      +3.4 | 73,697.862 | 64,493.702 |
| 106 | RID-GNSS-104  | 99.761524 | 14.598111 |  +40.5 |      +3.3 | 73,812.386 | 64,528.996 |
| 107 | RID-GNSS-105  | 99.766607 | 14.578337 |  +39.6 |      +3.6 | 74,361.400 | 62,341.374 |
| 108 | RID-GNSS-106  | 99.767685 | 14.579784 |  +40.3 |      +3.5 | 74,477.507 | 62,501.578 |
| 109 | RID-GNSS-107  | 99.768618 | 14.548066 |  +50.3 |      +2.0 | 74,580.070 | 58,992.174 |
| 110 | RID-GNSS-108  | 99.768535 | 14.549066 |  +49.3 |      +2.2 | 74,571.059 | 59,102.791 |
| 111 | RID-GNSS-109  | 99.723825 | 14.525301 |  +63.6 |      -1.5 | 69,753.606 | 56,470.993 |
| 112 | RID-GNSS-110  | 99.725094 | 14.523757 |  +64.5 |      -1.6 | 69,890.401 | 56,300.204 |
| 113 | RID-GNSS-111  | 99.762054 | 14.523128 |  +59.8 |      +0.3 | 73,874.209 | 56,232.493 |
| 114 | RID-GNSS-112  | 99.762451 | 14.525825 |  +58.8 |      +0.4 | 73,916.886 | 56,530.923 |
| 115 | RID-GNSS-113  | 99.761134 | 14.509647 |  +67.2 |      -0.9 | 73,775.877 | 54,740.862 |
| 116 | RID-GNSS-114  | 99.763471 | 14.508717 |  +67.3 |      -0.9 | 74,027.837 | 54,638.174 |
| 117 | RID-GNSS-115  | 99.789258 | 14.505498 |  +73.1 |      -0.7 | 76,807.786 | 54,283.677 |
| 118 | RID-GNSS-116  | 99.789027 | 14.503554 |  +74.1 |      -0.9 | 76,782.978 | 54,068.606 |
| 119 | RID-GNSS-117  | 99.709484 | 14.514695 |  +67.6 |      -2.5 | 68,208.231 | 55,296.951 |
| 120 | RID-GNSS-118  | 99.709453 | 14.516464 |  +67.7 |      -2.5 | 68,204.884 | 55,492.712 |
| 121 | RID-GNSS-119  | 99.691728 | 14.498472 |  +73.6 |      -3.8 | 66,294.752 | 53,501.362 |
| 122 | RID-GNSS-120  | 99.690255 | 14.499156 |  +76.0 |      -4.2 | 66,135.965 | 53,576.996 |
| 123 | RID-GNSS-121  | 99.723370 | 14.492606 |  +74.1 |      -3.2 | 69,705.946 | 52,853.491 |
| 124 | RID-GNSS-122  | 99.725132 | 14.492070 |  +71.9 |      -2.8 | 69,895.890 | 52,794.232 |
| 125 | RID-GNSS-123  | 99.774028 | 14.489141 |  +70.6 |      -1.0 | 75,167.161 | 52,472.791 |
| 126 | RID-GNSS-124  | 99.774295 | 14.487756 |  +71.2 |      -1.1 | 75,196.014 | 52,319.625 |
| 127 | RID-GNSS-125  | 99.693078 | 14.477380 |  +69.3 |      -3.1 | 66,440.954 | 51,167.741 |
| 128 | RID-GNSS-126  | 99.692867 | 14.475732 |  +68.1 |      -2.9 | 66,418.274 | 50,985.362 |
| 129 | RID-GNSS-127  | 99.714963 | 14.426613 |  +41.5 |      +1.7 | 68,802.319 | 45,551.410 |
| 130 | RID-GNSS-128  | 99.716913 | 14.426584 |  +41.0 |      +1.8 | 69,012.553 | 45,548.264 |
| 131 | RID-GNSS-129  | 99.750341 | 14.466167 |  +56.5 |      +0.4 | 72,614.899 | 49,929.473 |
| 132 | RID-GNSS-130  | 99.749555 | 14.464015 |  +55.9 |      +0.4 | 72,530.312 | 49,691.325 |
| 133 | RID-GNSS-131  | 99.752495 | 14.408349 |  +39.0 |      +3.2 | 72,850.514 | 43,532.331 |
| 134 | RID-GNSS-132  | 99.752163 | 14.406261 |  +38.5 |      +3.2 | 72,814.823 | 43,301.373 |
| 135 | RID-GNSS-133  | 99.674629 | 14.436152 |  +68.2 |      -3.2 | 64,452.765 | 46,605.647 |
| 136 | RID-GNSS-134  | 99.674669 | 14.434049 |  +72.3 |      -3.8 | 64,457.156 | 46,372.936 |
| 137 | RID-GNSS-135  | 99.719665 | 14.405527 |  +35.5 |      +2.8 | 69,310.206 | 43,218.573 |
| 138 | RID-GNSS-136  | 99.718165 | 14.406385 |  +36.2 |      +2.6 | 69,148.402 | 43,313.465 |
| 139 | RID-GNSS-137  | 99.667409 | 14.409900 |  +60.9 |      -2.1 | 63,674.761 | 43,700.874 |
| 140 | RID-GNSS-138  | 99.669360 | 14.409892 |  +63.7 |      -2.5 | 63,885.140 | 43,700.088 |
| 141 | RID-GNSS-139  | 99.684880 | 14.380455 |  +41.5 |      +1.2 | 65,559.550 | 40,443.333 |
| 142 | RID-GNSS-140  | 99.686511 | 14.379477 |  +40.4 |      +1.4 | 65,735.481 | 40,335.221 |
| 143 | RID-GNSS-141  | 99.652371 | 14.415809 |  +71.5 |      -3.9 | 62,052.965 | 44,354.536 |
| 144 | RID-GNSS-142  | 99.652317 | 14.414112 |  +69.6 |      -3.6 | 62,047.112 | 44,166.778 |
| 145 | RID-GNSS-143  | 99.669022 | 14.381103 |  +48.6 |      -0.2 | 63,849.147 | 40,514.810 |
| 146 | RID-GNSS-144  | 99.670737 | 14.380113 |  +47.4 |      +0.0 | 64,034.125 | 40,405.206 |
| 147 | RID-GNSS-145  | 99.644939 | 14.397687 |  +69.5 |      -3.6 | 61,251.577 | 42,349.406 |
| 148 | RID-GNSS-146  | 99.645064 | 14.396117 |  +67.6 |      -3.3 | 61,265.159 | 42,175.671 |
| 149 | RID-GNSS-147  | 99.681027 | 14.360542 |  +40.8 |      +1.2 | 65,144.450 | 38,240.041 |
| 150 | RID-GNSS-148  | 99.682250 | 14.360305 |  +41.6 |      +1.1 | 65,276.353 | 38,213.916 |
| 151 | RID-GNSS-149  | 99.621412 | 14.373073 |  +75.1 |      -4.5 | 58,714.154 | 39,626.021 |
| 152 | RID-GNSS-150  | 99.621746 | 14.372025 |  +77.0 |      -4.8 | 58,750.212 | 39,510.069 |
| 153 | RID-GNSS-151  | 99.638069 | 14.358826 |  +70.3 |      -3.8 | 60,510.785 | 38,049.680 |
| 154 | RID-GNSS-152  | 99.638122 | 14.357459 |  +69.3 |      -3.6 | 60,516.567 | 37,898.392 |
| 155 | RID-GNSS-153  | 99.595585 | 14.362570 |  +80.2 |      -5.1 | 55,928.337 | 38,464.254 |
| 156 | RID-GNSS-154  | 99.594318 | 14.362464 |  +80.6 |      -5.1 | 55,791.684 | 38,452.511 |
| 157 | RID-GNSS-155  | 99.635049 | 14.339189 |  +56.4 |      -1.6 | 60,185.066 | 35,876.939 |
| 158 | RID-GNSS-156  | 99.637192 | 14.338243 |  +57.0 |      -1.7 | 60,416.294 | 35,772.297 |
| 159 | RID-GNSS-157  | 99.609680 | 14.326041 |  +66.4 |      -3.1 | 57,448.251 | 34,422.369 |
| 160 | RID-GNSS-158  | 99.611089 | 14.325989 |  +65.3 |      -2.9 | 57,600.227 | 34,416.574 |
| 161 | RID-GNSS-159  | 99.582010 | 14.337042 |  +88.2 |      -6.2 | 54,463.532 | 35,640.011 |
| 162 | RID-GNSS-160  | 99.583411 | 14.338557 |  +85.7 |      -5.8 | 54,614.648 | 35,807.611 |
| 163 | RID-GNSS-161  | 99.600523 | 14.311107 |  +85.1 |      -5.9 | 56,460.148 | 32,770.211 |
| 164 | RID-GNSS-162  | 99.598708 | 14.311169 |  +87.7 |      -6.3 | 56,264.334 | 32,777.047 |
| 165 | RID-GNSS-163  | 99.665222 | 14.299128 |  +39.3 |      +1.2 | 63,440.544 | 31,444.747 |
| 166 | RID-GNSS-164  | 99.664828 | 14.297822 |  +38.8 |      +1.3 | 63,398.111 | 31,300.294 |
| 167 | RID-GNSS-165  | 99.608313 | 14.294164 |  +76.5 |      -4.6 | 57,300.377 | 30,895.439 |
| 168 | RID-GNSS-166  | 99.608019 | 14.293116 |  +78.7 |      -5.0 | 57,268.658 | 30,779.517 |
| 169 | RID-GNSS-167  | 99.635288 | 14.294844 |  +51.2 |      -0.8 | 60,210.938 | 30,970.559 |
| 170 | RID-GNSS-168  | 99.636712 | 14.294827 |  +50.1 |      -0.6 | 60,364.575 | 30,968.625 |
| 171 | RID-GNSS-169  | 99.539682 | 14.290166 |  +95.4 |      -6.4 | 49,895.228 | 30,455.010 |
| 172 | RID-GNSS-170  | 99.539020 | 14.288440 |  +94.7 |      -6.3 | 49,823.746 | 30,264.061 |
| 173 | RID-GNSS-171  | 99.543372 | 14.254574 |  +74.9 |      -3.3 | 50,291.837 | 26,516.823 |
| 174 | RID-GNSS-172  | 99.544038 | 14.252100 |  +77.5 |      -3.7 | 50,363.640 | 26,243.053 |
| 175 | RID-GNSS-173  | 99.521159 | 14.263374 |  +77.8 |      -3.1 | 47,895.268 | 27,491.558 |
| 176 | RID-GNSS-174  | 99.522098 | 14.261514 |  +75.5 |      -2.8 | 47,996.501 | 27,285.757 |
| 177 | RID-GNSS-175  | 99.520699 | 14.213118 |  +71.0 |      -2.0 | 47,842.916 | 21,931.158 |
| 178 | RID-GNSS-176  | 99.521716 | 14.213376 |  +71.6 |      -2.2 | 47,952.713 | 21,959.674 |





