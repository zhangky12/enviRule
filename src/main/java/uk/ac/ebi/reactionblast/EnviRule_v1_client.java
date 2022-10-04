package uk.ac.ebi.reactionblast;
import java.net.*;
import java.io.*;

public class EnviRule_v1_client {

    private Socket socket            = null;
    private DataInputStream  input   = null;
    private DataOutputStream out     = null;
    private DataInputStream res = null;

    // constructor to put ip address and port
    public EnviRule_v1_client(String address, int port)
    {
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
        }
        catch(UnknownHostException u)
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
        while (!line.equals("Over"))
        {
            try
            {
                if(line.equals("start")) start=true;

                line = input.readLine();

                if(line.equals("Over")) continue;

                out.writeUTF(line);
                out.flush();

                if(start){
//                    String respond = res.readUTF();
                    BufferedReader bs = new BufferedReader(
                            new InputStreamReader(socket.getInputStream()));
                    System.out.println("Start already");
                    System.out.println(bs.readLine());
                }

            }
            catch(IOException i)
            {
                System.out.println(i);
            }
        }

        // close the connection
        try
        {
            input.close();
            out.close();
            socket.close();
        }
        catch(IOException i)
        {
            System.out.println(i);
        }
    }

    public static void main(String args[])
    {
        EnviRule_v1_client client = new EnviRule_v1_client("127.0.0.1", 5000);
    }
}
