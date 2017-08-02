package test;

public class TimeUtils {
    private static double startTime;
    private static double endTime;

    public static void tic() {
        startTime = System.currentTimeMillis();
    }

    public static void toc() {
        endTime = System.currentTimeMillis();
        System.out.println((endTime - startTime) / 1000 + " 秒");
    }

    public static void toc(StringBuilder str) {
        endTime = System.currentTimeMillis();
        str.append((endTime - startTime) / 1000 + " 秒");
    }
}
