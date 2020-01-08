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

import java.awt.*;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.datatransfer.*;
import java.util.*;
import javax.swing.JDialog;
import java.io.*;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import java.awt.Dimension;

class ParamTable extends JPanel{
	public JTable table;
	public ParamTableModel tableModel;
	int nMax;
	private JScrollPane paramScrollPane;
	
	
	public ParamTable() {
		super(new GridLayout(1,1));
		tableModel = new ParamTableModel(1,4);
		table = new JTable(tableModel);
		table.setPreferredScrollableViewportSize(new Dimension(157, 97));
		table.setDragEnabled(true);
		JTableHeader tHeader;
		tHeader = table.getTableHeader();
		Color saumon = new Color( 255, 204, 153 );
		tHeader.setBackground(saumon);		
		Dimension dim = new Dimension( 157, 97);
		paramScrollPane = new JScrollPane(table);
		paramScrollPane.setSize(dim);
		nMax = 0;
		add(paramScrollPane);
		setBorder(BorderFactory.createEmptyBorder(13,0,0,0));	
				
	}
	public String serializeData(){
		String fileName = "para.txt";
		FileWriter fwData;
		StringBuffer sBuf  = new StringBuffer();
		
		try
		{
			int imax = 4;
			System.out.println(fileName);			
			
			fwData = new FileWriter(fileName);
			sBuf.append(nMax +  "\r\n");
			for( int i = 0; i < nMax ; i++)

				{
					sBuf.append(tableModel.data[i][0] + " " + tableModel.data[i][1]);
					if( ((Boolean)tableModel.data[i][2]).booleanValue() )
						sBuf.append(" 1");
					else 
						sBuf.append(" 0");
					if( ((Boolean)tableModel.data[i][3]).booleanValue() )
						sBuf.append(" 1");
					else 
						sBuf.append(" 0");
					sBuf.append("\r\n");
					
				}

				fwData.write(sBuf.toString(), 0 , sBuf.toString().length());
				fwData.close();
		}
		catch (IOException ee)
		{
			System.out.println(ee.toString());
		}
		return sBuf.toString();

	}
	
	public void importData(){
		nMax = 0;
		String fileName = "para.txt";	
		try
			{
				FileReader s;
				int i = 0;
				s = new FileReader(fileName);
				BufferedReader b = new BufferedReader(s);
				String content;
				content = b.readLine();			
				int nbLine = Integer.parseInt(content);
				tableModel = null; 
				tableModel = new ParamTableModel(nbLine,4);
				table.setModel(tableModel);
				double d;

							
				for( int k = 0  ; k  < nMax ;k ++){
					tableModel.data[k][0] = new String();
					tableModel.data[k][1] = new String();
					tableModel.data[k][2] = new Boolean(false);
					tableModel.data[k][3] = new Boolean(false);
				}

				while( content != null )
				{
					content = b.readLine();
					if ( content == null)
						break;
					StringTokenizer st = new StringTokenizer(content);
					if (st.countTokens() < 2) 
						break;		
					tableModel.data[i][0] = st.nextToken();
					tableModel.data[i][1] = st.nextToken();
					d = Double.parseDouble(st.nextToken() );
					if ( d < 0.5  ) 
						tableModel.data[i][2] = new Boolean(false);
					else 
						tableModel.data[i][2] = new Boolean(true);
					d = Double.parseDouble(st.nextToken() );
					if ( d < 0.5  ) 
						tableModel.data[i][3] = new Boolean(false);
					else 
						tableModel.data[i][3] = new Boolean(true);


					i++;
				}
				nMax = i;

				b.close();
				s.close();

			}
			catch ( IOException eIO	)
			{
				System.out.println(eIO.toString());
			}

	}
	public void setParamTableEnabled(boolean state)
	{
		paramScrollPane.setEnabled(state);
		table.setEnabled(state);
	}

	public void reinitParamTable(){
	tableModel = null; 
	tableModel = new ParamTableModel(1,4);
	table.setModel(tableModel);	
	}

	
}




class ParamTableModel extends AbstractTableModel
{
	int rows, cols;
	    
    private String[] columnNames = {"N","Freq","Effect"};		
    public Object[][] data;
	public ParamTableModel( int row, int col){
		rows = row;
		cols = col;
		
		int i, j ;
		data = new Object[rows][cols];
		for( i =0 ; i < rows  ; i++){
			for (j = 0 ; j < cols -2 ; j++ )
			{
				data[i][j] = new String();
			}
			for( ; j < cols ; j++) 
				data[i][j] = new Boolean(false);
		}


	}

	public void resetData()
	{
		for( int i =0 ; i < rows ; i++)
			for (int j  = 0 ; j < cols ; j++ )
			{
				data[i][j] = new String();
			}
	}
	public void initData(Object d[][], int a, int b)
	{
		for( int i = 0 ; i < a ; i++)
			for (int j  = 0 ; j < b ; j++ )
			{
				data[i][j] = (String)d[i][j];
			}
	}

	public int getColumnCount() {
		return columnNames.length;
	}

	public int getRowCount() {
		return data.length;
	}

	public String getColumnName(int col) {
		return columnNames[col];
	}

	public Object getValueAt(int row, int col) {
		if( data == null)
			return null;
		return data[row][col];
	}


	public Class getColumnClass(int c) {
		return getValueAt(0, c).getClass();
	}


	public boolean isCellEditable(int row, int col) {
		if( col == 0 || col == 1 )
			return false;
		return true;
	}
	

	public void setValueAt(Object value, int row, int col) {
		data[row][col] = value;
		fireTableCellUpdated(row, col);

	}

};
