function plot_data(data)
    s = size(data,2);
    figure
    subplot(s,1,1)

    plot(data(:,1),data(:,2),'b')
    hold on;
    plot(data(:,1),data(:,end),'r--')
    for i = 3:1:s
        subplot(s,1,i-1)
        plot(data(:,1),data(:,i))
        
    end
    
end
