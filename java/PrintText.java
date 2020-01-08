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
import java.awt.font.*;
import java.awt.geom.*;
import java.awt.print.*;
import java.text.*;


public class PrintText implements Printable {


	private AttributedString mStyledText;
	private String mText;
	private String sPrint;
	public PrintText(String s)
	{
		sPrint = s;
		mStyledText = new AttributedString(s);
	}


	public void runPrint() {


		PrinterJob printerJob = PrinterJob.getPrinterJob();


		Book book = new Book();
		book.append(new PrintText(sPrint), new PageFormat());


		printerJob.setPageable(book);



		boolean doPrint = printerJob.printDialog();
		if (doPrint) {


		try {

			printerJob.print();

		} catch (PrinterException exception) {

			System.err.println("Printing error: " + exception);

		}

		}

		}

		public int print(Graphics g, PageFormat format, int pageIndex) {


		Graphics2D g2d = (Graphics2D) g;


		g2d.translate(format.getImageableX(), format.getImageableY());

		g2d.setPaint(Color.black);

		Point2D.Float pen = new Point2D.Float();
		AttributedCharacterIterator charIterator = mStyledText.getIterator();
		LineBreakMeasurer measurer = new LineBreakMeasurer(charIterator, g2d.getFontRenderContext());
		float wrappingWidth = (float) format.getImageableWidth();


		while (measurer.getPosition() < charIterator.getEndIndex()) {

		TextLayout layout = measurer.nextLayout(wrappingWidth);
		pen.y += layout.getAscent();
		float dx = layout.isLeftToRight()? 0 : (wrappingWidth - layout.getAdvance());

		layout.draw(g2d, pen.x + dx, pen.y);
		pen.y += layout.getDescent() + layout.getLeading();

		}
		return Printable.PAGE_EXISTS;

	}
}
