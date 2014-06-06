using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.IO;
namespace SFproject
{
    /// <summary>
    /// MainWindow.xaml \u7684\u4EA4\u4E92\u903B\u8F91
    /// </summary>
    public partial class MainWindow : Window
    {
        [DllImport("Project1.dll",
                EntryPoint = "getsearch",
                CharSet = CharSet.Ansi,
                  CallingConvention = CallingConvention.StdCall)]
        public static extern int getsearch(int[] a, [MarshalAs(UnmanagedType.LPStr)]string f);
        public static int[] result = new int[450];
        public static int[] lenseq = new int[450];
        public static  double MAX = 10000000;
        public static int ENV = 10;
        public static string filename = "";
        private MediaPlayer player = new MediaPlayer();
        private static double dtw(double[] a, int lena, double[] b, int lenb)
        {
            //求lena和lenb的最小公倍数
           
           int min_len = 100;  //lena和lenb的最小公倍数
           min_len = Math.Min(lena, lenb);
           ENV = Convert.ToInt32(Math.Ceiling((min_len) * 0.01));
           double[] a1 = new double[min_len];
           double[] b1 = new double[min_len];
           if (lena < lenb)
           {

               for (int i = 1; i <= min_len; i++)
               {
                   double num = 0;
                   int j = lenb * (i - 1) / min_len + 1;
                   int max = lenb * i / min_len;
                   int delta = max - j+1;
                   for (; j <= max; j++)
                   {
                       num += b[j - 1];
                   }
                   b1[i - 1] = num / delta;
                   //cout << now[i-1] << endl;
               }
               for (int i = 0; i < min_len; i++)
                   a1[i] = a[i];

           }
           else
           {
               for (int i = 1; i <= min_len; i++)
               {
                   double num = 0;
                   int j = lena * (i - 1) / min_len + 1;
                   int max = lena * i / min_len;
                   int delta = max - j + 1;
                   for (; j <= max; j++)
                   {
                       num += a[j - 1];
                   }
                   a1[i - 1] = num / delta;
                   //cout << now[i-1] << endl;
               }
               for (int i = 0; i < min_len; i++)
                   b1[i] = b[i];
           }
          

          /*  int a_temp = lena, b_temp = lenb;
            if (a_temp > b_temp)       //a>b 互换a b值
            {
                a_temp = a_temp + b_temp;
                b_temp = a_temp - b_temp;
                a_temp = a_temp - b_temp;
            }
            for (min_len = b_temp; ; min_len += b_temp)
            {
                if (min_len % a_temp == 0 && min_len % b_temp == 0)//满足最小公倍数条件
                    break;
            }
            //	cout << min_len << endl;
            //升采样   到相同长度
            double[] a1 = new double[min_len];
            int wa = min_len / lena;
            //	cout << wa << endl;
            for (int i = 1; i <= min_len; i++)
            {
                a1[i - 1] = a[(int)(Math.Ceiling((double)i / (double)wa)) - 1];
                //cout << ceil(float(i)/float(wa)) << " ";
            }
            //cout << endl;
            double[] b1 = new double[min_len];
            int wb = min_len / lenb;
            //	cout << wb << endl;
            for (int i = 1; i <= min_len; i++)
            {
                b1[i - 1] = b[(int)(Math.Ceiling((double)i / (double)wb)) - 1];
            }*/

            /*	for(int i = 0;i < min_len;i++)
                    cout << a1[i] <<" ";
                cout << endl;
                for(int i = 0;i < min_len;i++)
                    cout << b1[i] << " ";
                cout << endl;*/
            //动规求LDTW
            //动态创建二维数组
            double[,] memory = new double[min_len, min_len];

            //初始化
            memory[min_len - 1, min_len - 1] = Math.Pow((double)(a1[min_len - 1] - b1[min_len - 1]), 2);

            for (int i = min_len - 2; i >= 0; i--)
            {
                if (min_len - 1 - i <= ENV)
                {
                    memory[min_len - 1, i] = Math.Pow((double)(a1[min_len - 1] - b1[i]), 2) + memory[min_len - 1, i + 1];
                    memory[i, min_len - 1] = Math.Pow((double)(a1[i] - b1[min_len - 1]), 2) + memory[i + 1, min_len - 1];

                }
                else
                {
                    memory[min_len - 1, i] = MAX;
                    memory[i, min_len - 1] = MAX;

                }
            }
            //计算
            for (int i = min_len - 2; i >= 0; i--)
            {
                double min = memory[i, i + 1];
                if (memory[i + 1, i + 1] < min)
                    min = memory[i + 1, i + 1];
                if (memory[i + 1, i] < min)
                    min = memory[i + 1, i];
                memory[i, i] = Math.Pow(a1[i] - b1[i], 2) + min;
                for (int j = i - 1; j >= 0; j--)
                {
                    if (i - j <= ENV)
                    {
                        min = memory[i + 1, j + 1];

                        if (memory[i + 1, j] < min)
                            min = memory[i + 1, j];
                        if (memory[i, j + 1] < min)
                            min = memory[i, j + 1];
                        memory[i, j] = Math.Pow(a1[i] - b1[j], 2) + min;

                        min = memory[j + 1, i + 1];
                        if (memory[j + 1, i] < min)
                            min = memory[j + 1, i];
                        if (memory[j, i + 1] < min)
                            min = memory[j, i + 1];
                        memory[j, i] = Math.Pow(a1[j] - b1[i], 2) + min;

                    }
                    else
                    {
                        memory[i, j] = MAX;
                        memory[j, i] = MAX;
                    }
                }
            }
            double result = memory[0, 0];
            
            return result/min_len;
        }
        private static Dictionary<int, string> musicnamedic = new Dictionary<int, string>();//编号与歌曲名
        private static Dictionary<int, double[]> seqdic = new Dictionary<int, double[]>();//编号与原始序列
        private static Process p;
        
        public MainWindow()
        {
            InitializeComponent();
            Stopbut.IsEnabled = false;
            Pausebut.IsEnabled = false;
           StreamReader sr=new StreamReader("h:\\projectdata\\seq_full.txt");
           probar.Visibility = Visibility.Hidden;
           string tmp;
           int id = 0;
          
           while ((tmp = sr.ReadLine())!=null)//歌曲名
           {
               musicnamedic.Add(id, tmp);
               tmp = sr.ReadLine();
               lenseq[id] = Convert.ToInt32(tmp);//长度
               tmp = sr.ReadLine();
               string[] sarray = tmp.Split();
               double[] darray=new double[sarray.Length];//获取序列
               int count=0;
               double sum = 0;
               foreach (string s in sarray)
               {
                   darray[count] = Convert.ToDouble(s);
                   sum += darray[count];
                   count++;
               }
               sum = sum / lenseq[id];
               for (int i = 0; i < count; i++)
               {
                   darray[i] -= sum;
               }
               seqdic.Add(id, darray);
               id++;
              
           }
           sr.Close();
        }

        private void onwave(object sender, RoutedEventArgs e)
        {
            System.Windows.Forms.OpenFileDialog openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            openFileDialog1.InitialDirectory ="D:\\Documents\\Music";
            openFileDialog1.Filter = "music files (*.mp3,*.wav,*.wma)|*.mp3;*.wav;*.wma";
            openFileDialog1.FilterIndex = 0;
            openFileDialog1.RestoreDirectory = true;
            if (openFileDialog1.ShowDialog() == System.Windows.Forms.DialogResult.OK)
            {
                probar.Visibility = Visibility.Visible;
                Filename.Text = openFileDialog1.FileName;
                filename = "h:\\projectdata\\temp.txt";
                string command = "Extract.exe " + openFileDialog1.FileName + " > " + "h:\\projectdata\\temp.txt";
                p = new Process();
                p.Exited += p_Exited;
                p.EnableRaisingEvents = true;
                p.StartInfo.FileName = "cmd.exe";
                p.StartInfo.UseShellExecute = false;
                p.StartInfo.RedirectStandardInput = true;
                p.StartInfo.RedirectStandardOutput = true;
                p.StartInfo.RedirectStandardError = true;
                p.StartInfo.CreateNoWindow = true;
                p.Start();
                p.StandardInput.WriteLine(command);//"Extract.exe h:\\我的歌声里1.wma > h:\\projectdata\\temp.txt"
                p.StandardInput.WriteLine("exit");
              
                
            }
           
        }
        private void p_Exited(object sender, EventArgs e)
        {
            this.probar.Dispatcher.Invoke(
                new Action(
                    delegate
                    {
                        probar.Visibility = Visibility.Hidden;
                    })
                    );
        }
        private void onsearch(object sender, RoutedEventArgs e)
        {
            p.WaitForExit();
            media.Stop();
            int preid=0;
            /*if(!predict.Text.Equals(""))
             preid = Convert.ToInt32(predict.Text);*/
            Stopbut.IsEnabled = false;
            Pausebut.IsEnabled = false;
            //new string("h:\\projectdata\\下雨天.txt")
            
            int hitnum = getsearch(result, filename);
            StreamReader sr = new StreamReader(filename);
            string tmp;
            tmp = sr.ReadLine();
             
            int length = Convert.ToInt32(tmp);
            tmp = sr.ReadLine();
            string[] sarray = tmp.Split();
           
            double[] darray = new double[sarray.Length];
            int count = 0;
            double sum = 0;
            foreach (string s in sarray)
            {
                darray[count] = Convert.ToDouble(s);
                sum += darray[count];
                count++;
            }
            sum = sum / length;
            for (int i = 0; i < count; i++)
            {
                darray[i] -= sum;
            }
            dis_id[] dislist = new dis_id[hitnum];
            for (int i = 0; i < hitnum; i++)
            {
                double[] testseq = seqdic[result[i]];
                double distance=dtw(darray,length,testseq,lenseq[result[i]]);
                
                dislist[i].distance = distance;
                dislist[i].id=result[i];
                
            }
            Array.Sort(dislist);
           // Array.Sort(dislist, new test());
           // hitnum=Math.Min(hitnum,20);
            Viewerpan.Children.Clear();
            for (int i = 0; i < hitnum; i++)
            {
                Button bt = new Button();
                bt.Tag = (int)0;
                bt.Click += playasong;
                bt.Content = musicnamedic[dislist[i].id];
               /* if (preid != 0)
                {
                    if (dislist[i].id == preid)
                        predict.Text = "at " + i.ToString();
                }*/
                
                Viewerpan.Children.Add(bt);

            }
            sr.Close();
          /*  media.Stop();
            
            media.Source = new Uri("E:\\BaiduMUSIC\\Songs\\年度之歌 - 谢安琪.mp3", UriKind.Absolute);
            media.Position = TimeSpan.Zero;
            media.Play();*/
            //player.Source = new Uri("E:\\BaiduMUSIC\\Songs\\年度之歌 - 谢安琪.mp3", UriKind.Absolute);
         
            
        }

        private void onmediafailed(object sender, ExceptionRoutedEventArgs e)
        {
            System.Windows.MessageBox.Show(e.ErrorException.Message);
        }
        private void onpause(object sender, RoutedEventArgs e)
        {
            Button bt = sender as Button;
            if (bt != null)
            {
                if ((string)bt.Content == "Pause")
                {
                    media.Pause();
                    bt.Content = "Go on";
                }
                else
                {
                    media.Play();
                    bt.Content = "Pause";
                }
            }
        }
        private void onstop(object sender, RoutedEventArgs e)
        {
            media.Stop();
        }
        private void playasong(object sender, RoutedEventArgs e)
        {

            media.Stop();
            Pausebut.Content = "Pause";
            Button bt = sender as Button;
            if (bt != null)
            {
                Stopbut.IsEnabled = true;
                Pausebut.IsEnabled = true;
                string name = (string)bt.Content;
                string path = "E:\\BaiduMUSIC\\Songs\\songlib\\";
                path += name;
                path += ".mp3";
                media.Source = new Uri(path, UriKind.Absolute);
                media.Position = TimeSpan.Zero;
                media.Play();


            }
        }
    }
}
public struct dis_id : IComparable<dis_id>
{
    public double distance;
    public int id;
    public int CompareTo(dis_id x)
    {
        if (distance > x.distance)
            return 1;
        else return -1;

    }

    
}
public class test : IComparer<dis_id>
{
    public int Compare(dis_id x,dis_id y)
    {
        if (x.distance > y.distance)
            return 1;
        else return -1;

    }
}
