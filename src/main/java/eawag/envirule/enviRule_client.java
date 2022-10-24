package eawag.envirule;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.Socket;
import java.net.UnknownHostException;

public class enviRule_client {

    private Socket socket            = null;
    private DataInputStream input   = null;
    private DataOutputStream out     = null;
    private DataInputStream res = null;

    public enviRule_client(String address, int port) {

        // establish a connection
        try
        {
            socket = new Socket(address, port);
            System.out.println("Connected");

            // takes input from terminal
            input  = new DataInputStream(System.in);

            // sends output to the socket
            out    = new DataOutputStream(socket.getOutputStream());

            // receives respond from the socket
//            res = new DataInputStream(socket.getInputStream());
        }catch(UnknownHostException u)
        {
            System.out.println(u);
        }
        catch(IOException i)
        {
            System.out.println(i);
        }

        // string to read message from input
        String line = "";

        // keep reading until "Over" is input
        boolean start = false;



    }
}
