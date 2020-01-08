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

 public class ColumnHeaderListener extends MouseAdapter {
        public void mouseClicked(MouseEvent evt) {
            JTable table = ((JTableHeader)evt.getSource()).getTable();
            TableColumnModel colModel = table.getColumnModel();
    
            // The index of the column whose header was clicked
            int vColIndex = colModel.getColumnIndexAtX(evt.getX());
            int mColIndex = table.convertColumnIndexToModel(vColIndex);
   
	    
			TableColumn tm = table.getColumnModel().getColumn(0);
			  tm.setCellRenderer(new ColorColumnRenderer(Color.lightGray, Color.blue));

            // Return if not clicked on any column header
            if (vColIndex == -1) {
                return;
            }
    
            // Determine if mouse was clicked between column heads
            Rectangle headerRect = table.getTableHeader().getHeaderRect(vColIndex);
            if (vColIndex == 0) {
				
				
				
                headerRect.width -= 3;    
            } else {
				
				
				
                headerRect.grow(-3, 0);   
            }
            if (!headerRect.contains(evt.getX(), evt.getY())) {
                int vLeftColIndex = vColIndex;
                if (evt.getX() < headerRect.x) {
                    vLeftColIndex--;
                }
            }
        }

    }
	class ColorColumnRenderer extends DefaultTableCellRenderer 
	{
	   Color bkgndColor, fgndColor;

		public ColorColumnRenderer(Color bkgnd, Color foregnd) {
			super(); 
			bkgndColor = bkgnd;
			fgndColor = foregnd;
		}

		public Component getTableCellRendererComponent
						(JTable table, Object value, boolean isSelected,
						 boolean hasFocus, int row, int column) 
		{
			Component cell = super.getTableCellRendererComponent( table, value, isSelected,
																  hasFocus, row, column);

			cell.setBackground( bkgndColor );
			cell.setForeground( fgndColor );
			return cell;
		}
	}
