/*
 * The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
 * Copyright (C) 1998-2000 University of South Carolina
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 * 
 * The current version can be found at <A
 * HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>
 * 
 * Bug reports and comments should be directed to H. Philip Crotwell,
 * crotwell@seis.sc.edu or Tom Owens, owens@seis.sc.edu
 * 
 */
package edu.sc.seis.TauP;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.Socket;

/**
 * example TauP client. Connects to the TauPDaemon and calculates times.
 * 
 * @version 1.1.3 Wed Jul 18 15:00:35 GMT 2001
 * 
 * 
 * 
 * @author H. Philip Crotwell
 * 
 */
public class TauPClient {

    public static void main(String[] args) {
        DataInputStream fromServer;
        DataOutputStream toServer;
        try {
            Socket timeSocket = new Socket("localhost", 6371, true);
            fromServer = new DataInputStream(timeSocket.getInputStream());
            toServer = new DataOutputStream(timeSocket.getOutputStream());
            toServer.writeBytes("MODEL prem\n");
            toServer.writeBytes("PHASES P,S\n");
            toServer.writeBytes("DEPTH 100.0\n");
            toServer.writeBytes("DISTANCE 60.0\n");
            String outString;
            while((outString = fromServer.readLine()) != null) {
                System.out.println(outString);
            }
        } catch(IOException e) {}
    }
}
