/*
 * This file is part of THESIAS
 * Copyright (C) 2004-2020 David-Alexandre Trégouët, Valérie Garelle
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

public class thesiaslib
{
	private String FileName;
	private int Nbloci;
	private int []Idloci;
	private	int Ldmatrix;
	private int Msdata;
	private int R2;
	private int Chxt;
	private int Num0;
	private int Idtime;
	private int Idoffset;
	private int Ajust;
	private int [] Numajust;
	private int Maxvarfic;
	
	private int Offset;
	private int Xlnk;
	private int Numsx;
	
		
	public thesiaslib()
	{
		FileName = null;
		Nbloci   = 0;
		Idloci	 = null;
		Ldmatrix = 0;
		Msdata	 = 0;
		R2	 = 0;
		Chxt	 = 0;
		Num0	 = 0;
		Idtime	 = 0;
		Offset   = 0;
		Idoffset  = 0;
		Ajust	 = 0;
		Numajust = null;
		Maxvarfic = 0;
		Xlnk	  = 0;
		Numsx     = 0;
	}

	public thesiaslib(String sfilename,
			  int nbloci,
			  int [] idloci,
			  int ldmatrix,
			  int msdata,
			  int pR2,
			  int chxt,
			  int num0,
			  int idtime,
			  int idoffset,
			  int ajust,
			  int [] numajust,
			  int maxvarfic )
	{
		FileName = sfilename;
		Nbloci   = nbloci;
		Idloci	 = idloci;
		Ldmatrix = ldmatrix;
		Msdata	 = msdata;
		R2	 = pR2;
		Chxt	 = chxt;
		Num0	 = num0;
		Idtime	 = idtime;
		Idoffset = idoffset;
		Ajust	 = ajust;
		Numajust = numajust;
		Maxvarfic = maxvarfic;
	}

	public void SetThesias(String sfilename,
			       int nbloci,
			       int [] idloci,
			       int ldmatrix,
			       int msdata,
			       int pR2,
			       int chxt,
			       int num0,
			       int idtime,
			       int offset,
			       int idoffset,
			       int ajust,
			       int [] numajust,
			       int maxvarfic,
			       int xlnk,
			       int numsx)
	{
		FileName = sfilename;
		Nbloci   = nbloci;
		Idloci	 = idloci;
		Ldmatrix = ldmatrix;
		Msdata	 = msdata;
		R2	 = pR2;
		Chxt	 = chxt;
		Num0	 = num0;
		Idtime	 = idtime;
		Offset   = offset;
		Idoffset = idoffset;
		Ajust	 = ajust;
		Numajust = numajust;
		Maxvarfic = maxvarfic;
		Xlnk	  = xlnk;
		Numsx   = numsx;
	}
	
	public void SetThesias(String sfilename,
			       int nbloci,
			       int [] idloci,
			       int ldmatrix,
			       int msdata,
			       int pR2,
			       int chxt,
			       int num0,
			       int idtime,
			       int idoffset,
			       int ajust,
			       int [] numajust,
			       int maxvarfic)
	{
		FileName = sfilename;
		Nbloci   = nbloci;
		Idloci	 = idloci;
		Ldmatrix = ldmatrix;
		Msdata	 = msdata;
		R2	 = pR2;
		Chxt	 = chxt;
		Num0	 = num0;
		Idtime	 = idtime;
		Idoffset = idoffset;
		Ajust	 = ajust;
		Numajust = numajust;
		Maxvarfic = maxvarfic;
	}

	public void run()
	{	

		System.out.println("  FileName : " + FileName);
		System.out.println("  Maxvarfic: " + Maxvarfic);
		System.out.println("  Nbloci   : " + Nbloci);
		for ( int j = 0; j < Nbloci; j++ )
		{
			System.out.println("	Idloci ["+ j +"]   : " + Idloci[j]);
		}
		System.out.println("  Ldmatrix : " + Ldmatrix);
		System.out.println("  Msdata   : " + Msdata);
		System.out.println("  R2       : " + R2);
		System.out.println("  Chxt     : " + Chxt);
		System.out.println("  Num0     : " + Num0);
		System.out.println("  Idtime   : " + Idtime);
		System.out.println("  Offset   : " + Offset);
		System.out.println("  Idoffset : " + Idoffset);
		System.out.println("  Ajust    : " + Ajust);
		for ( int j = 0; j < Ajust; j++ )
		{
			System.out.println("	Numajust ["+ j +"] : " + Numajust[j]);
		}
		System.out.println("  Xlnk     : " + Xlnk);
		System.out.println("  Numsx    : " + Numsx);

		int i = -222;
		try
		{
			i = thesiasRun( FileName,
				        Maxvarfic,
				        Nbloci,
					Idloci,
					Ldmatrix,
					Msdata,
					R2,
					Chxt,
					Num0,
					Idtime,
					Offset,
					Idoffset,
					Ajust,
					Numajust,
					Xlnk,
					Numsx);
		}
		catch(Exception e)
		{
			System.out.println(e.toString());
		}
	}

	public native int  thesiasRun(String sFileName,
				      int maxvarfic,
				      int nbloci,
				      int [] idloci,
				      int ldmatrix,
				      int msdata,
				      int R2,
				      int chxt,
				      int num0,
				      int idtime,
				      int offset,
				      int idoffset,
				      int ajust,
				      int [] numajust,
				      int xlnk,
				      int numsx);

	static
	{System.loadLibrary("thesiaslib");}

}
