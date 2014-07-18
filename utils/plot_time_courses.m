function plot_time_courses(X, plot_envs, channel_string, target_source_idx)


if not(exist('plot_envs', 'var'))
    plot_envs = 1;
end
if not(exist('channel_string', 'var'))
    channel_string = 'x';
end
if not(exist('target_source_idx', 'var'))
    target_source_idx = [];
end

N = size(X,2);
stds = std(X);

Env = abs(hilbert(X));

hold on
box on
y_offset = [0 cumsum(5*stds(1:end-1))];
y_labels = cell(1,N);
for k=1:N
    plot(X(:,k)+y_offset(k));
    if plot_envs
        if any(ismember(k, target_source_idx))
            env_color = 'r';
        else
            env_color = 0.6*[1,1,1];
        end
        plot(Env(:,k) + y_offset(k), 'color', env_color)
    end
    y_labels{k} = sprintf('%s_%d',channel_string,k);
end
ylim([-4, y_offset(end)+4])
set(gca, 'ytick', y_offset)
set(gca, 'ytickLabel', y_labels);